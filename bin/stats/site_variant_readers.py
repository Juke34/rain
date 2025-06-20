import numpy as np
from numpy.typing import NDArray
from typing import TextIO, Optional
from utils import SiteVariantData, NUC_STR_TO_IND, EDIT_TYPES, NONEDIT_TYPES
from abc import ABC, abstractmethod

REDITOOLS_FIELDS = ["Seqid", "Position", "Reference", "Strand", "Coverage", "MeanQ", "Frequencies"]
REDITOOLS_FIELD_INDEX = {field: i for i, field in  enumerate(REDITOOLS_FIELDS)}

class RNAVariantReader(ABC):
    """Abstract class defining the API for readers"""

    @abstractmethod
    def __init__(self, file_handle: TextIO) -> None:
        pass

    @abstractmethod
    def read(self) -> Optional[SiteVariantData]:
        pass

    @abstractmethod
    def close(self) -> None:
        pass


class TestReader(RNAVariantReader):
    def __init__(self, strand: int, edit: str) -> None:
        """Create a TestReader that returns SiteVariantData objects with only one read for one type of edition.
        Arguments:
        - strand: Strand of the simulated features
        - edit: Kind of edit (or non-edit) to simulate. Examples: AA, AC, AG, AT, CC, CA. 
        """
        self.position: int = 0
        self.strand:int = strand
        assert edit in EDIT_TYPES + NONEDIT_TYPES
        self.reference: int = NUC_STR_TO_IND.get(edit[0], 4)
        self.edited: str = edit[1]
        self.frequencies: NDArray[np.int32] = np.zeros(5, dtype=np.int32)
        self.frequencies[NUC_STR_TO_IND.get(self.edited, 4)] = 1

        return None

    def read(self) -> Optional[SiteVariantData]:
        data = SiteVariantData(
            seqid="test",
            position=self.position,
            reference=self.reference,
            strand=self.strand,
            coverage=1,
            mean_quality=30.0,
            frequencies=self.frequencies
        )
        self.position += 1

        return data
        
    def close(self) -> None:
        pass


class ReditoolsXReader(RNAVariantReader):
    header_strings = (
        "Region",
        "Position",
        "Reference",
        "Strand",
        "Coverage-q30",
        "MeanQ",
        "BaseCount[A,C,G,T]",
        "AllSubs",
        "Frequency",
        "gCoverage-q30",
        "gMeanQ",
        "gBaseCount[A,C,G,T]",
        "gAllSubs",
        "gFrequency",
    )

    def __init__(self, file_handle: TextIO) -> None:
        self.file_handle: TextIO = file_handle

        line = self.file_handle.readline()
        if line == "":
            raise Exception("Empty file")

        self._get_parts(line)
        self._validate_header()

        self.file_position: int = 1
        self.genomic_position: int = 0

        return None

    def _get_parts(self, line: str) -> None:
        self.parts: list[str] = [s.strip() for s in line.split("\t")]

        return None

    def _validate_header(self) -> None:
        l1 = len(self.header_strings)
        l2 = len(self.parts)

        if l1 != l2:
            raise Exception(
                f"The header of the input file contains {l2} fields instead of {l1} fields."
            )

        for a, b in zip(self.header_strings, self.parts):
            if a != b:
                raise Exception(
                    f'Unexpected field name "{b}" in the header of the input file (expected "{a}")'
                )

        return None
    
    def parse_strand(self):
        pass

    def _parse_parts(self) -> SiteVariantData:
        strand = self.parse_strand()
            
        reference_nuc_str: str = self.parts[REDITOOLS_FIELD_INDEX["Reference"]]
        
        return SiteVariantData(
            seqid=self.parts[REDITOOLS_FIELD_INDEX["Seqid"]],
            position=int(self.parts[REDITOOLS_FIELD_INDEX["Position"]]) - 1,    # Convert Reditools 1-based index to Python's 0-based index
            reference=NUC_STR_TO_IND.get(reference_nuc_str, 4),
            strand=strand,
            coverage=int(self.parts[REDITOOLS_FIELD_INDEX["Coverage"]]),
            mean_quality=float(self.parts[REDITOOLS_FIELD_INDEX["MeanQ"]]),
            frequencies=np.int32(self.parts[REDITOOLS_FIELD_INDEX["Frequencies"]][1:-1].split(",") + [0])
        )

    def read(self) -> Optional[SiteVariantData]:
        """Read the data of the next variant site"""
        line = self.file_handle.readline()

        if line == "":
            # End of file reached
            return None

        self._get_parts(line)

        return self._parse_parts()
    
    def close(self) -> None:
        """Close the file"""
        self.file_handle.close()
    

class Reditools2Reader(ReditoolsXReader):
    def parse_strand(self) -> int:
        strand = int(self.parts[REDITOOLS_FIELD_INDEX["Strand"]])
        match strand:
            case 0:
                return -1
            case 1:
                return 1
            case 2:
                return 0
            case _:
                raise Exception(f"Invalid strand value: {strand}")

class Reditools3Reader(ReditoolsXReader):
    def parse_strand(self) -> int:
        strand_str = self.parts[REDITOOLS_FIELD_INDEX["Strand"]]
        match strand_str:
            case "-":
                return -1
            case "+":
                return 1
            case "*":
                return 0
            case _:
                raise Exception(f"Invalid strand value: {strand_str}")

