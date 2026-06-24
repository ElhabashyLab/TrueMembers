import logging
import os
from io import StringIO

from Bio import SeqIO

from True_Members.scripts.utils import extract_uid_from_recordID

logger = logging.getLogger(__name__)


class UniprotLocalDB:
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.index_path = fasta_path + ".idx"
        self._index: dict[str, int] | None = None


    def _build_index(self):
        logger.info("Building index for %s ...", self.fasta_path)
        index: dict[str, int] = {}
        with open(self.fasta_path, "r") as f:
            offset = f.tell()
            line = f.readline()
            while line:
                if line.startswith(">"):
                    uid = extract_uid_from_recordID(line)
                    if uid:
                        index[uid] = offset
                offset = f.tell()
                line = f.readline()

        tmp_path = self.index_path + f".tmp.{os.getpid()}"
        with open(tmp_path, "w") as f:
            for uid, off in index.items():
                f.write(f"{uid}\t{off}\n")
        os.rename(tmp_path, self.index_path)
        logger.info("Index built: %d entries at %s", len(index), self.index_path)
        return index


    def _load_index(self) -> dict[str, int]:
        if not os.path.exists(self.index_path):
            return self._build_index()

        index: dict[str, int] = {}
        with open(self.index_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                uid, off_str = line.split("\t", 1)
                index[uid] = int(off_str)
        logger.info("Loaded index: %d entries from %s", len(index), self.index_path)
        return index


    def ensure_index(self):
        if self._index is None:
            self._index = self._load_index()


    def get_sequence(self, accession: str):
        self.ensure_index()
        offset = self._index.get(accession)
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
