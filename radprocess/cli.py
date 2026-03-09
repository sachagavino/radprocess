from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.align import Align
from rich import box
from rich.table import Table
import sys
from importlib.metadata import version

import socket
import getpass

from radprocess.interface.interface import launch_interface

console = Console()

VERSION = "0.1.0"
DOC_URL = "https://radprocess.readthedocs.io/en/latest/"


def print_banner():

    title = Text("RadProcess", style="bold orange1")
    VERSION = version("radprocess")

    body = Text()
    body.append("Welcome to RadProcess\n", style="grey82")
    body.append(
        "A post-processing pipeline for multi-fluid RAMSES simulations.\n\n",
        style="grey70",
    )

    body.append("Documentation: ", style="bold")
    body.append(f"{DOC_URL}\n\n", style="underline bright_blue")

    body.append(
        "A Gradio interface will open below. If a public share link is available, you can access it remotely.",
        style="dark_blue",
    )

    panel = Panel(
        Align.center(body),
        title=title,
        subtitle=f"v{VERSION}",
        box=box.DOUBLE,
        padding=(1, 4),
        border_style="bright_blue",
    )

    console.print(panel)


def print_status_table():

    table = Table(title="Environment", box=box.SIMPLE)

    table.add_column("Item", style="bold")
    table.add_column("Value", style="cyan")

    table.add_row("radprocess version", version("radprocess"))
    table.add_row("Python", sys.version.split()[0])
    table.add_row("User", getpass.getuser())
    table.add_row("Host", socket.gethostname())

    console.print(table)


def print_ssh_panel(port=7860):

    user = getpass.getuser()
    host = socket.gethostname()

    text = Text()

    text.append("Gradio share link could not be created.\n\n", style="bold yellow")

    text.append(
        "This usually happens on secured servers or HPC clusters where "
        "outbound SSH tunnels are blocked.\n\n",
        style="grey70",
    )

    text.append(
        "To access the RadProcess interface from your laptop run:\n\n",
        style="bold",
    )

    text.append(
        f"ssh -L {port}:localhost:{port} {user}@{host}\n\n",
        style="bold cyan",
    )

    text.append(
        "Then open the following address in your browser:\n\n",
        style="bold",
    )

    text.append(
        f"http://localhost:{port}",
        style="bold green",
    )

    panel = Panel(
        Align.center(text),
        title="[bold orange1]RadProcess Connection Instructions[/bold orange1]",
        box=box.DOUBLE,
        border_style="bright_blue",
        padding=(1, 4),
    )

    console.print(panel)


def main():

    if sys.version_info < (3, 11):
        console.print("[bold red]Error:[/bold red] radprocess requires Python ≥ 3.11")
        sys.exit(1)

    print_banner()

    print_status_table()

    console.print()

    port = 7860

    try:

        with console.status(
            "[bold green]Launching RadProcess interface...[/bold green]",
            spinner="dots",
        ):
            launch_interface(share=True)

    except Exception:

        console.print(
            "\n[yellow]Could not create Gradio share link. "
            "Falling back to local server.[/yellow]\n"
        )

        print_ssh_panel(port)

        with console.status(
            "[bold green]Starting local RadProcess server...[/bold green]",
            spinner="dots",
        ):
            launch_interface(
                share=False,
                server_name="0.0.0.0",
                server_port=port,
            )