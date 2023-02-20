from pathlib import Path

from livereload import Server, shell

if __name__ == "__main__":
    server = Server()
    server.watch(
        str(Path(__file__).resolve().parent / "*.py") + "/*",
        shell("make.bat html"),
        delay=1,
    )
    server.watch(
        str(Path(__file__).resolve().parent / "*.rst") + "/*",
        shell("make.bat html"),
        delay=1,
    )
    server.watch(
        str(Path(__file__).resolve().parent / "*.md") + "/*",
        shell("make.bat html"),
        delay=1,
    )
    server.watch(
        str(Path(__file__).resolve().parent / "*.ipynb") + "/*",
        shell("make.bat html"),
        delay=1,
    )
    server.serve(root="build/html")