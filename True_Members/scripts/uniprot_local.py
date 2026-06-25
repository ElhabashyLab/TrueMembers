import logging
import os
import subprocess
from io import StringIO

from Bio import SeqIO

from True_Members.scripts.utils import extract_uid_from_recordID

logger = logging.getLogger(__name__)


class UniprotLocalDB:
    """Provides fast local lookup of UniProt sequences by accession.

    Indexes a FASTA file by accession→byte offset and uses binary search for
    O(log n) lookup, avoiding repeated full-file scans or network calls.
    """

    def __init__(self, fasta_path: str):
        """Initialize with a path to the UniProt FASTA file.

        Args:
            fasta_path: Path to a FASTA file containing UniProt sequences.
        """
        self.fasta_path = fasta_path
        self.index_path = fasta_path + ".idx"

    def _build_index(self):
        """Build a sorted tab-separated index file (accession → byte offset).

        Scans the FASTA file header lines, extracts UniProt accessions, and
        writes a sorted index so binary search can locate entries quickly.
        Uses the C locale to guarantee sort order matches Python's string
        comparison used in :meth:`_binary_search_offset`.
        """
        logger.info("Building index for %s ...", self.fasta_path)
        tmp_path = self.index_path + f".tmp.{os.getpid()}"
        count = 0
        with open(self.fasta_path, "r") as fa, open(tmp_path, "w") as out:
            offset = fa.tell()
            line = fa.readline()
            while line:
                if line.startswith(">"):
                    uid = extract_uid_from_recordID(line)
                    if uid:
                        out.write(f"{uid}\t{offset}\n")
                        count += 1
                offset = fa.tell()
                line = fa.readline()

        logger.info("Scanned %d entries, sorting index ...", count)
        subprocess.run(
            ["sort", "-o", self.index_path, tmp_path],
            check=True,
            env={**os.environ, "LC_ALL": "C"},
        )
        os.remove(tmp_path)
        logger.info("Index built: %d entries → %s", count, self.index_path)

    def ensure_index(self):
        """Build the index if it does not already exist.

        Safe to call before every lookup; only builds on first access.
        """
        if not os.path.exists(self.index_path):
            self._build_index()

    def _binary_search_offset(self, accession: str) -> int | None:
        """Binary-search the sorted index file for *accession*'s byte offset.

        Args:
            accession: UniProt accession string (e.g. ``P12345``).

        Returns:
            Byte offset into the FASTA file where the record header begins,
            or ``None`` if *accession* is not found.
        """
        with open(self.index_path, "r") as f:
            lo, hi = 0, os.path.getsize(self.index_path)
            while lo < hi:
                mid = (lo + hi) // 2
                f.seek(mid)
                while mid > 0:
                    f.seek(mid - 1)
                    if f.read(1) == "\n":
                        f.seek(mid)
                        break
                    mid -= 1
                    f.seek(mid)
                pos = f.tell()
                line = f.readline()
                if not line:
                    break
                try:
                    acc, off_str = line.strip().split("\t", 1)
                except ValueError:
                    lo = f.tell()
                    continue
                if acc == accession:
                    return int(off_str)
                elif acc < accession:
                    lo = f.tell()
                else:
                    hi = pos
        return None

    def get_sequence(self, accession: str):
        """Return the Biopython ``SeqRecord`` for *accession*, or ``None``.

        Builds the index on first call (if missing), then binary-searches for
        the byte offset and reads the record directly from the FASTA file.

        Args:
            accession: UniProt accession string.

        Returns:
            A :class:`Bio.SeqRecord.SeqRecord` or ``None`` if the accession
            is not present in the FASTA file.
        """
        self.ensure_index()
        offset = self._binary_search_offset(accession)
        if offset is None:
            return None

        with open(self.fasta_path, "r") as f:
            f.seek(offset)
            lines = []
            line = f.readline()
            while line:
                if line.startswith(">") and lines:
                    break
                lines.append(line)
                line = f.readline()

        return next(SeqIO.parse(StringIO("".join(lines)), "fasta"), None)
