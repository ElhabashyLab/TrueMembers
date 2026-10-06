import sys

from truemembers import __version__


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in ("--version", "-V"):
        print(f"truemembers {__version__}")
        return 0
    print("truemembers: CLI skeleton, commands not implemented yet", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
