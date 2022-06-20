import pysam
from scripts.tools import sequence as jseq
from collections import defaultdict
from tqdm import tqdm
import pickle
import click

@click.group()
def main():
    pass

def parseBam(ADD_SEQ_BAM, OUT_PICKLE, BY_PRIMER):

    addSeqBam = pysam.AlignmentFile(ADD_SEQ_BAM, "r")

    seqTransformer = jseq()
    bamDict = defaultdict()

    for x in tqdm(addSeqBam, "parse bam file", total=addSeqBam.mapped, ascii=True):
        if not BY_PRIMER:
            readESSeq = x.get_tag("ES")
            readFSSeq = x.get_tag("FS")
            if x.is_reverse:
                readStrand = 1
            else:
                readStrand = 0
            bamDict[x.qname]=(
                [
                    readESSeq,
                    seqTransformer.reverseComplement(readESSeq),
                    readFSSeq,
                    seqTransformer.reverseComplement(readFSSeq),
                    readStrand,
                ]
            )
        else:
            readPSSeq = x.get_tag("PS")
            if x.is_reverse:
                readStrand = 1
            else:
                readStrand = 0
            bamDict[x.qname]=(
                [
                    readPSSeq,
                    seqTransformer.reverseComplement(readPSSeq),
                    readStrand,
                ]
            )

    with open(OUT_PICKLE, 'wb') as handle:
        pickle.dump(bamDict, handle)


@main.command("parseBam")
@click.option("-b", "ADD_SEQ_BAM", help="bam added unmapped seq tag")
@click.option("-o", "OUT_FEATHER", help="output feather format")
@click.option(
    "--by-primer",
    "BY_PRIMER",
    is_flag=True,
    help="Extraction of region between primers and aligned sequences")

def _parseBam( ADD_SEQ_BAM, OUT_FEATHER, BY_PRIMER ):
    parseBam( ADD_SEQ_BAM, OUT_FEATHER, BY_PRIMER )

main()
