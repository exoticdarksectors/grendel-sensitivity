"""gzip a file the way the tracked senscalc_2501/*.gz were written (macOS `gzip -n -9`, which is zlib-based):
zlib deflate level 9, memLevel 8, header mtime 0, XFL 2, OS 3 (Unix). GNU gzip -n -9 uses its own deflate and gives
different bytes for the same content; Python's gzip module writes OS 255.

    python3 -I gzip_zlib.py <in> <out.gz>
"""
import struct
import sys
import zlib


def main() -> int:
    data = open(sys.argv[1], "rb").read()
    c = zlib.compressobj(9, zlib.DEFLATED, -15, 8, zlib.Z_DEFAULT_STRATEGY)
    body = c.compress(data) + c.flush()
    head = b"\x1f\x8b\x08\x00" + b"\x00\x00\x00\x00" + b"\x02\x03"
    tail = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data) & 0xFFFFFFFF)
    open(sys.argv[2], "wb").write(head + body + tail)
    return 0


if __name__ == "__main__":
    sys.exit(main())
