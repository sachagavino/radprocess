"""
Tests for radprocess core components.

These tests verify the fundamental building blocks of the pipeline:
physical constants, configuration dataclasses, grid construction,
OcTree structure, dust size distributions, and RAMSES file parsing.
Designed to run without external simulation data (RAMSES, RADMC-3D).

Note: pymses is mocked in conftest.py so tests can run in CI
without pymses installed.
"""

import os
import tempfile

import numpy as np
import pytest

from radprocess.constants.constants import (
    c, kbol, M_sun, L_sun, R_sun, T_sun, au2cm, au2m, pc2cm,
    amu, proton_mass_g, Na, sigma, dust_to_gas,
)
from radprocess.utils.config import (
    Directories, AmrSource, Sim, Radmc3dConfig, PolarisConfig,
    ConfigParams, fancy_repr, fmt_value,
)
from radprocess.utils.ramsesinfo import SinkInfo, html_table_from_columns
from radprocess.pipeline.Grid import Grid
from radprocess.pipeline.OcTree import OcTree, CellOct
from radprocess.dustproperties.SizeDistrib import SizeDistrib
from radprocess.ramses.read import get_snapshot_number, hydro_file_descriptor


# ============================================================
# Physical constants
# ============================================================

class TestConstants:
    """Verify physical constants are consistent and reasonable."""

    def test_speed_of_light(self):
        assert pytest.approx(c, rel=1e-3) == 2.998e10

    def test_solar_mass(self):
        assert pytest.approx(M_sun, rel=1e-2) == 1.989e33

    def test_solar_luminosity(self):
        assert pytest.approx(L_sun, rel=1e-1) == 3.8525e33

    def test_stefan_boltzmann(self):
        assert pytest.approx(sigma, rel=1e-3) == 5.670e-5

    def test_au_to_cm(self):
        assert pytest.approx(au2cm, rel=1e-3) == 1.496e13

    def test_au_to_m(self):
        assert pytest.approx(au2m, rel=1e-3) == au2cm * 1e-2

    def test_parsec_to_cm(self):
        assert pytest.approx(pc2cm, rel=1e-3) == 3.086e18

    def test_avogadro_times_amu(self):
        """Na * amu should be ~1 g/mol."""
        assert pytest.approx(Na * amu, rel=1e-2) == 1.0

    def test_boltzmann_cgs(self):
        assert pytest.approx(kbol, rel=1e-2) == 1.381e-16

    def test_dust_to_gas_default(self):
        assert dust_to_gas == 0.01

    def test_solar_temperature(self):
        assert pytest.approx(T_sun, rel=1e-2) == 5780


# ============================================================
# Configuration dataclasses
# ============================================================

class TestDirectories:

    def test_defaults(self):
        d = Directories()
        assert d.ramses_output == "ramses_outputs/"
        assert d.pipeline_output == "pipeline_outputs/"

    def test_custom_values(self):
        d = Directories(ramses_output="/data/run1/", pipeline_output="/data/out/")
        assert d.ramses_output == "/data/run1/"

    def test_type_enforcement(self):
        with pytest.raises(TypeError):
            Directories(ramses_output=123)


class TestAmrSource:

    def test_defaults(self):
        a = AmrSource()
        assert a.rho is True
        assert a.vel is False
        assert a.dustratios is False

    def test_enable_fields(self):
        a = AmrSource(vel=True, bl=True, br=True)
        assert a.vel is True
        assert a.bl is True

    def test_type_enforcement(self):
        with pytest.raises(TypeError):
            AmrSource(rho="yes")


class TestSim:

    def test_defaults(self):
        s = Sim()
        assert s.size_hole_au == 0.0
        assert s.dtogas == 0.01
        assert s.facc == 0.1
        assert s.use_multi_grain is True

    def test_custom_dtogas(self):
        s = Sim(dtogas=0.001)
        assert s.dtogas == 0.001


class TestRadmc3dConfig:

    def test_defaults(self):
        r = Radmc3dConfig()
        assert r.nphot == 1_000_000
        assert r.setthreads == 8
        assert r.wave_min < r.wave_max
        assert r.n_wavelengths == 200

    def test_wavelength_range(self):
        r = Radmc3dConfig()
        assert r.wave_min < 1.0
        assert r.wave_max > 1000


class TestPolarisConfig:

    def test_defaults(self):
        p = PolarisConfig()
        assert p.dust_size_min < p.dust_size_max
        assert p.dust_size_powerlaw == -3.5
        assert p.mass_fraction == 1


class TestConfigParams:

    def test_defaults(self):
        cfg = ConfigParams()
        assert isinstance(cfg.dir, Directories)
        assert isinstance(cfg.amrsource, AmrSource)
        assert isinstance(cfg.sim, Sim)
        assert isinstance(cfg.radmc3d, Radmc3dConfig)
        assert isinstance(cfg.polaris, PolarisConfig)
        assert cfg.nb_dust == 0

    def test_repr_does_not_crash(self):
        cfg = ConfigParams()
        text = repr(cfg)
        assert len(text) > 0
        assert "PARAMS" in text

    def test_html_repr_does_not_crash(self):
        cfg = ConfigParams()
        html = cfg._repr_html_()
        assert "<details" in html


# ============================================================
# Utility functions
# ============================================================

class TestFmtValue:

    def test_boolean_true(self):
        assert fmt_value(True) == "True"

    def test_boolean_false(self):
        assert fmt_value(False) == "False"

    def test_large_number(self):
        result = fmt_value(1.5e10)
        assert "e" in result or "E" in result

    def test_string(self):
        result = fmt_value("hello")
        assert "hello" in result


class TestFancyRepr:

    def test_produces_output(self):
        d = Directories()
        text = fancy_repr(d)
        assert "Directories" in text
        assert "ramses_output" in text

    def test_non_dataclass_falls_back(self):
        result = fancy_repr("not a dataclass")
        assert result == "'not a dataclass'"


# ============================================================
# SinkInfo
# ============================================================

class TestSinkInfo:

    def test_creation(self):
        info = SinkInfo(
            columns=["id", "mass", "x", "y", "z"],
            rows=[{"id": "1", "mass": "0.5", "x": "0.0", "y": "0.0", "z": "0.0"}],
            num_sinks=1,
            data=np.array([[1, 0.5, 0.0, 0.0, 0.0]]),
        )
        assert info.num_sinks == 1
        assert info.data.shape == (1, 5)

    def test_repr(self):
        info = SinkInfo(
            columns=["id", "mass"],
            rows=[{"id": "1", "mass": "0.5"}],
            num_sinks=1,
            data=np.array([[1, 0.5]]),
        )
        assert "SinkTable" in repr(info)

    def test_html_repr(self):
        info = SinkInfo(
            columns=["id", "mass"],
            rows=[{"id": "1", "mass": "0.5"}],
            num_sinks=1,
            data=np.array([[1, 0.5]]),
        )
        assert "Sink File" in info._repr_html_()


class TestHtmlTableFromColumns:

    def test_dict_rows(self):
        html = html_table_from_columns(
            columns=["a", "b"], rows=[{"a": 1, "b": 2}], title="Test",
        )
        assert "<table" in html

    def test_list_rows(self):
        html = html_table_from_columns(
            columns=["x", "y"], rows=[[10, 20]], title="Table",
        )
        assert "<table" in html


# ============================================================
# Grid
# ============================================================

class TestGrid:

    def test_init_empty(self):
        g = Grid()
        assert g.density == []
        assert g.stars == []
        assert g.amr_grid == []

    def test_add_density(self):
        g = Grid()
        g.add_density(np.ones((10, 10, 10)))
        assert len(g.density) == 1
        assert g.density[0].shape == (10, 10, 10)

    def test_add_star(self):
        g = Grid()
        g.add_star({"mass": 1.0, "temp": 5780})
        assert len(g.stars) == 1

    def test_add_temperature(self):
        g = Grid()
        g.add_temperature(np.full((5, 5, 5), 300.0))
        assert len(g.temperature) == 1

    def test_wavelength_grid_linear(self):
        g = Grid()
        g.set_wavelength_grid(lmin=0.1, lmax=1000, nlam=100)
        assert len(g.lam) == 100
        assert g.lam[0] == pytest.approx(0.1)
        assert g.lam[-1] == pytest.approx(1000)

    def test_wavelength_grid_log(self):
        g = Grid()
        g.set_wavelength_grid(lmin=0.1, lmax=1000, nlam=100, log=True)
        dl_short = g.lam[1] - g.lam[0]
        dl_long = g.lam[-1] - g.lam[-2]
        assert dl_short < dl_long

    def test_mono_wavelength_grid(self):
        g = Grid()
        g.set_mcmonowavelength_grid(lmin=0.5, lmax=2.0, nlam=50, log=True)
        assert len(g.monolam) == 50

    def test_add_multiple_densities(self):
        g = Grid()
        for i in range(3):
            g.add_density(np.ones((5, 5, 5)) * i)
        assert len(g.density) == 3
        assert g.density[2].mean() == pytest.approx(2.0)

    def test_add_amr_grid(self):
        g = Grid()
        g.add_amr_grid({"type": "octree", "levels": 5})
        assert len(g.amr_grid) == 1


# ============================================================
# OcTree
# ============================================================

class TestCellOct:

    def test_creation(self):
        cell = CellOct(0.0, 0.0, 0.0, 1.0, 0)
        assert cell.x_min == 0.0
        assert cell.length == 1.0
        assert cell.level == 0
        assert cell.isleaf == 0
        assert cell.data == []
        assert cell.branches == []


class TestOcTree:

    def test_creation(self):
        tree = OcTree(0.0, 0.0, 0.0, 1.0)
        assert tree.root is not None
        assert tree.root.length == 1.0
        assert tree.cell_counter == 0

    def test_init_cell_boundaries(self):
        """Splitting a cell should create 8 children with half the parent length."""
        tree = OcTree(0.0, 0.0, 0.0, 1.0)
        tree.initCellBoundaries(tree.root, 1)
        assert len(tree.root.branches) == 8
        for child in tree.root.branches:
            assert child.length == pytest.approx(0.5)
            assert child.level == 1

    def test_child_positions_cover_parent(self):
        """Children should tile the parent cell — 8 unique positions."""
        tree = OcTree(0.0, 0.0, 0.0, 2.0)
        tree.initCellBoundaries(tree.root, 1)
        positions = set()
        for child in tree.root.branches:
            positions.add((child.x_min, child.y_min, child.z_min))
        assert len(positions) == 8

    def test_insert_leaf(self):
        """Inserting a cell at the correct level should make it a leaf."""
        tree = OcTree(0.0, 0.0, 0.0, 1.0)
        leaf = CellOct(0.0, 0.0, 0.0, 0.5, 1)
        leaf.data = [1.0, 2.0, 3.0]
        tree.initCellBoundaries(tree.root, 1)
        tree.insertInTree(tree.root, leaf, 0)
        inserted = tree.root.branches[0]
        assert inserted.isleaf == 1
        assert inserted.data == [1.0, 2.0, 3.0]

    def test_check_octree_leaf(self):
        tree = OcTree(0.0, 0.0, 0.0, 1.0)
        tree.root.isleaf = 1
        tree.root.data = [0.0]
        tree.nr_of_cells = 1
        assert tree.checkOcTree(tree.root) is True

    def test_check_octree_none(self):
        tree = OcTree(0.0, 0.0, 0.0, 1.0)
        assert tree.checkOcTree(None) is True

    def test_reset_counter(self):
        tree = OcTree(0.0, 0.0, 0.0, 1.0)
        tree.cell_counter = 42
        tree.reset_counter()
        assert tree.cell_counter == 0


# ============================================================
# SizeDistrib
# ============================================================

class TestSizeDistrib:

    def test_sizes_shape(self):
        sd = SizeDistrib(ndust=10, amin=0.005, amax=1000.0)
        sizes = sd.sizes()
        assert sizes.shape == (3, 10)

    def test_sizes_span_range(self):
        sd = SizeDistrib(ndust=5, amin=0.01, amax=100.0)
        sizes = sd.sizes()
        assert pytest.approx(sizes[0, 0], rel=1e-3) == 0.01
        assert pytest.approx(sizes[1, -1], rel=1e-3) == 100.0

    def test_average_between_edges(self):
        sd = SizeDistrib(ndust=5, amin=0.01, amax=100.0)
        sizes = sd.sizes()
        for i in range(5):
            assert sizes[0, i] < sizes[2, i] < sizes[1, i]

    def test_sizes_monotonically_increasing(self):
        sd = SizeDistrib(ndust=8, amin=0.005, amax=500.0)
        sizes = sd.sizes()
        for i in range(7):
            assert sizes[0, i] < sizes[0, i + 1]

    def test_log_spacing(self):
        """Bins should be log-spaced (constant ratio between consecutive edges)."""
        sd = SizeDistrib(ndust=4, amin=1.0, amax=10000.0)
        sizes = sd.sizes()
        ratios = sizes[0, 1:] / sizes[0, :-1]
        assert pytest.approx(ratios[0], rel=1e-6) == ratios[1]

    def test_single_bin_spans_full_range(self):
        sd = SizeDistrib(ndust=1, amin=0.1, amax=10.0)
        sizes = sd.sizes()
        assert sizes.shape == (3, 1)
        assert pytest.approx(sizes[0, 0], rel=1e-6) == 0.1
        assert pytest.approx(sizes[1, 0], rel=1e-6) == 10.0


# ============================================================
# RAMSES file parsing
# ============================================================

class TestGetSnapshotNumber:

    def test_standard_path(self):
        assert get_snapshot_number("output_00042") == "00042"

    def test_path_with_slash(self):
        assert get_snapshot_number("output_00100/") == "00100"

    def test_no_underscore(self):
        assert get_snapshot_number("nosnapshothere") is None

    def test_integer_input(self):
        assert get_snapshot_number(12345) is None

    def test_nested_path(self):
        assert get_snapshot_number("/data/sim/output_00001") == "00001"


class TestHydroFileDescriptor:

    def test_valid_descriptor(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            content = (
                "nvar = 5\n"
                "variable #  1: density\n"
                "variable #  2: velocity_x\n"
                "variable #  3: velocity_y\n"
                "variable #  4: velocity_z\n"
                "variable #  5: dust_ratio_1\n"
            )
            with open(os.path.join(tmpdir, "hydro_file_descriptor.txt"), "w") as f:
                f.write(content)
            nvar, variables, nb_dust = hydro_file_descriptor(tmpdir)
            assert nvar == 5
            assert len(variables) == 5
            assert variables[1] == "density"
            assert nb_dust == 1

    def test_multiple_dust_ratios(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            content = (
                "nvar = 8\n"
                "variable #  1: density\n"
                "variable #  2: velocity_x\n"
                "variable #  3: velocity_y\n"
                "variable #  4: velocity_z\n"
                "variable #  5: dust_ratio_1\n"
                "variable #  6: dust_ratio_2\n"
                "variable #  7: dust_ratio_3\n"
                "variable #  8: thermal_pressure\n"
            )
            with open(os.path.join(tmpdir, "hydro_file_descriptor.txt"), "w") as f:
                f.write(content)
            nvar, variables, nb_dust = hydro_file_descriptor(tmpdir)
            assert nvar == 8
            assert nb_dust == 3

    def test_missing_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(FileNotFoundError):
                hydro_file_descriptor(tmpdir)

    def test_no_dust_ratios(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            content = (
                "nvar = 4\n"
                "variable #  1: density\n"
                "variable #  2: velocity_x\n"
                "variable #  3: velocity_y\n"
                "variable #  4: velocity_z\n"
            )
            with open(os.path.join(tmpdir, "hydro_file_descriptor.txt"), "w") as f:
                f.write(content)
            nvar, variables, nb_dust = hydro_file_descriptor(tmpdir)
            assert nb_dust == 0
