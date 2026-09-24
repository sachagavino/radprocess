"""
radprocess.polaris.view_geometry
================================

Single source of truth for POLARIS detector viewing geometry.

The problem this solves
-----------------------
A POLARIS ``<detector_dust>`` map is oriented by two *rotation* axes
(``<axis1>``, ``<axis2>``) plus two angles (``rot_angle_1`` = theta,
``rot_angle_2`` = phi).  POLARIS builds the image frame by rotating the
identity basis by theta about axis1, then by phi about axis2
(see ``CDetector::setDetCoordSystem``).  The line of sight is the resulting
``ez``; the image plane is spanned by ``ex`` (horizontal) and ``ey`` (vertical).

The historical bug: ``axis1``/``axis2`` were also used, directly, as if they
were the image-plane basis -- both to hand-write the standard views and to
project the sink offset for detector centering.  They coincide with ``ex``/``ey``
only for ``xy`` and ``xz``; for ``yz`` they do not, so the line of sight and the
recentering shift were both wrong.

The fix: define every view by the thing we actually care about -- the image
axes ``right`` (=ex) and ``up`` (=ey) -- and *derive* the POLARIS rotation
parameters from that single frame.  The detector shift is always projected onto
``ex``/``ey``.  A round-trip assertion guarantees the parameters POLARIS receives
reproduce exactly the frame we intend, so the two can never silently disagree.
"""

from __future__ import annotations

import numpy as np

_ATOL = 1e-6


# ----------------------------------------------------------------------
#  Rotation primitives -- must match POLARIS Vector3D::rot exactly
#  (right-handed Rodrigues rotation about a unit axis).
# ----------------------------------------------------------------------
def _rodrigues(v, axis, angle_rad):
    v = np.asarray(v, dtype=float)
    n = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(n)
    if norm == 0.0:
        # POLARIS leaves the vector unchanged when the angle is 0; a zero
        # axis only ever reaches here paired with a zero angle.
        return v.copy()
    n = n / norm
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return v * c + np.cross(n, v) * s + n * np.dot(n, v) * (1.0 - c)


def frame_from_polaris(axis1, axis2, theta_deg, phi_deg):
    """Replicate ``CDetector::setDetCoordSystem``.

    Returns ``(ex, ey, ez)`` unit vectors: ``ex`` horizontal, ``ey`` vertical,
    ``ez`` the map normal (line of sight / observer direction).
    """
    t = np.radians(theta_deg)
    p = np.radians(phi_deg)
    out = []
    for v0 in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
        v = _rodrigues(v0, axis1, t)
        v = _rodrigues(v, axis2, p)
        out.append(v / np.linalg.norm(v))
    return out[0], out[1], out[2]


def polaris_from_frame(right, up):
    """Inverse: given the desired image axes, return POLARIS parameters.

    ``right`` -> ex (horizontal), ``up`` -> ey (vertical); the line of sight is
    ``ez = ex x up``.  Uses a single rotation (phi = 0), so ``axis1`` is the
    axis-angle axis of the frame and ``axis2`` is an unused placeholder.

    Returns ``(axis1, theta_deg, axis2, phi_deg, ex, ey, ez)``.
    """
    ex = np.asarray(right, dtype=float)
    ey = np.asarray(up, dtype=float)
    ex = ex / np.linalg.norm(ex)
    ey = ey / np.linalg.norm(ey)
    if abs(np.dot(ex, ey)) > 1e-9:
        raise ValueError(
            f"View axes must be orthogonal: right={right}, up={up} "
            f"(dot={np.dot(ex, ey):.3e})"
        )
    ez = np.cross(ex, ey)

    # Rotation matrix taking the identity basis to (ex, ey, ez).
    M = np.column_stack([ex, ey, ez])

    # Axis-angle decomposition of M (a proper rotation, det = +1).
    cos_t = np.clip((np.trace(M) - 1.0) / 2.0, -1.0, 1.0)
    theta = np.arccos(cos_t)

    if theta < 1e-12:
        axis = np.array([0.0, 0.0, 1.0])
        theta = 0.0
    elif abs(theta - np.pi) < 1e-9:
        # 180 deg: sin(theta) ~ 0, recover the axis from (M + I)/2.
        A = 0.5 * (M + np.eye(3))
        k = int(np.argmax(np.diag(A)))
        axis = A[:, k] / np.sqrt(max(A[k, k], 0.0))
    else:
        axis = np.array([
            M[2, 1] - M[1, 2],
            M[0, 2] - M[2, 0],
            M[1, 0] - M[0, 1],
        ]) / (2.0 * np.sin(theta))
    axis = axis / np.linalg.norm(axis)

    return axis, np.degrees(theta), np.array([0.0, 0.0, 1.0]), 0.0, ex, ey, ez


# ----------------------------------------------------------------------
#  Standard views -- defined by their IMAGE axes, not rotation axes.
#    right = ex (horizontal, increasing to the right)
#    up    = ey (vertical, increasing upward)
#    LOS   = ez = right x up
# ----------------------------------------------------------------------
STANDARD_VIEW_AXES = {
    "xy": {"right": [1, 0, 0], "up": [0, 1, 0], "plane_id": 1},  # LOS +z
    "xz": {"right": [1, 0, 0], "up": [0, 0, 1], "plane_id": 2},  # LOS -y
    "yz": {"right": [0, 1, 0], "up": [0, 0, 1], "plane_id": 3},  # LOS +x
}


def resolve_view_geometry(name, spec):
    """Return a validated geometry dict for one view.

    ``spec`` may be either
      * new form   : ``{"right": [...], "up": [...], "plane_id": int}``
      * legacy form: ``{"axis1": [...], "axis2": [...], "theta": deg,
                        "phi": deg, "plane_id": int}``

    The returned dict always contains ``axis1, axis2, theta, phi`` (what POLARIS
    receives) together with ``ex, ey, ez`` (the exact frame POLARIS builds from
    them, used for the detector shift).  A round-trip check guarantees the two
    descriptions agree.
    """
    plane_id = int(spec.get("plane_id", 1))

    if "right" in spec and "up" in spec:
        axis1, theta, axis2, phi, ex, ey, ez = polaris_from_frame(
            spec["right"], spec["up"]
        )
    elif all(k in spec for k in ("axis1", "axis2", "theta", "phi")):
        axis1 = np.asarray(spec["axis1"], dtype=float)
        axis2 = np.asarray(spec["axis2"], dtype=float)
        theta = float(spec["theta"])
        phi = float(spec["phi"])
        ex, ey, ez = frame_from_polaris(axis1, axis2, theta, phi)
    else:
        raise ValueError(
            f"View '{name}' must define either (right, up) or "
            f"(axis1, axis2, theta, phi); got keys {sorted(spec)}."
        )

    # Guarantee: the parameters POLARIS receives reproduce the frame we think.
    ex_c, ey_c, ez_c = frame_from_polaris(axis1, axis2, theta, phi)
    for got, want, lbl in ((ex_c, ex, "ex"), (ey_c, ey, "ey"), (ez_c, ez, "ez")):
        if not np.allclose(got, want, atol=_ATOL):
            raise RuntimeError(
                f"View '{name}': POLARIS rotation would build {lbl}={got}, "
                f"but the intended image basis is {want}. Refusing to render a "
                f"view whose geometry is internally inconsistent."
            )

    return {
        "plane_id": plane_id,
        "axis1": [float(x) for x in axis1],
        "axis2": [float(x) for x in axis2],
        "theta": float(theta),
        "phi": float(phi),
        "ex": [float(x) for x in ex],
        "ey": [float(x) for x in ey],
        "ez": [float(x) for x in ez],
    }


def detector_shift_2d(view_geom, shift_3d):
    """Project a 3-D offset (metres) onto the view's image plane (ex, ey)."""
    shift = np.asarray(shift_3d, dtype=float)
    ex = np.asarray(view_geom["ex"], dtype=float)
    ey = np.asarray(view_geom["ey"], dtype=float)
    return float(np.dot(shift, ex)), float(np.dot(shift, ey))


if __name__ == "__main__":
    # -------- self-test --------
    au = 1.495978707e11
    axis_name = ["x", "y", "z"]

    print("Standard views (derived POLARIS params + resulting LOS):")
    for name, spec in STANDARD_VIEW_AXES.items():
        g = resolve_view_geometry(name, spec)
        ez = np.round(g["ez"], 6)
        k = int(np.argmax(np.abs(ez)))
        los = ("+" if ez[k] >= 0 else "-") + axis_name[k]
        print(f"  {name}: axis1={np.round(g['axis1'],3)} theta={g['theta']:.3f} "
              f"phi={g['phi']:.1f} -> ex={np.round(g['ex'],3)} ey={np.round(g['ey'],3)} "
              f"ez(LOS)={ez} ({los})")

    # Legacy specs must still resolve and round-trip.
    legacy = {"axis1": [1, 0, 0], "axis2": [0, 0, 1], "theta": 90, "phi": 0}
    g = resolve_view_geometry("legacy_xz", legacy)
    assert np.allclose(g["ez"], [0, -1, 0], atol=_ATOL)

    # The real sink offset from the user's run; check yz recentering is fixed.
    off = np.array([-579.2, 421.7, -602.8]) * au
    g = resolve_view_geometry("yz", STANDARD_VIEW_AXES["yz"])
    dx, dy = detector_shift_2d(g, off)
    print("\nyz detector shift with the refactor:")
    print(f"  map_shift = ({dx:.6e}, {dy:.6e})")
    print(f"  = (offset_y, offset_z) = ({off[1]:.6e}, {off[2]:.6e})  -> centered on sink")
    assert np.allclose([dx, dy], [off[1], off[2]], atol=1.0)

    # A few random custom frames must round-trip exactly.
    rng = np.random.default_rng(0)
    for _ in range(2000):
        a = rng.normal(size=3)
        b = rng.normal(size=3)
        b = b - b.dot(a) / a.dot(a) * a          # make orthogonal
        if np.linalg.norm(a) < 1e-3 or np.linalg.norm(b) < 1e-3:
            continue
        g = resolve_view_geometry("rand", {"right": a, "up": b})
        ex, ey = np.array(a) / np.linalg.norm(a), np.array(b) / np.linalg.norm(b)
        assert np.allclose(g["ex"], ex, atol=_ATOL)
        assert np.allclose(g["ey"], ey, atol=_ATOL)
    print("\nAll round-trip and consistency checks passed.")