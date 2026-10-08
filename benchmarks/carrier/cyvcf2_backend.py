"""Benchmark-only reader adapter. Production modules are never edited."""
from types import SimpleNamespace
import numpy as np
from cyvcf2 import VCF
from extract_exact_carriers import ValidationError, text


class Reader:
    def __init__(self, path, index_filename):
        self.vcf = VCF(str(path), threads=1, lazy=False, strict_gt=False)
        try:
            self.vcf.set_index(str(index_filename))
            self.header = SimpleNamespace(contigs=self.vcf.seqnames, samples=self.vcf.samples)
            self.types = {h['ID']:dict(h.info()) for h in self.vcf.header_iter() if h.type == 'FORMAT'}
        except Exception:
            self.vcf.close()
            raise

    def fetch(self, contig, start, stop):
        for r in self.vcf(f'{contig}:{start+1}-{stop}'):
            yield Record(self, r)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.vcf.close()


class Record:
    def __init__(self, reader, record):
        self.reader, self.record = reader, record
        self.contig, self.pos, self.ref, self.alts = record.CHROM, record.POS, record.REF, record.ALT
        self.format = record.FORMAT


def format_value(array, sample, info):
    if array is None:
        return '.'
    value = array[sample]
    kind = info['Type']
    scalar = str(info['Number']) == '1'
    if kind in ('String', 'Character'):
        if isinstance(value, bytes):value=value.decode('utf-8')
        if isinstance(value,np.ndarray):value=value.tolist()
        if isinstance(value,list):value=','.join(map(str,value))
        return str(value)
    values = np.atleast_1d(value)
    result=[]
    for x in values:
        if kind == 'Integer':
            n=int(x)
            if n == -2147483647:break  # htslib vector-end padding, not a missing allele/depth
            result.append(None if n == -2147483648 else n)
        elif kind == 'Float':
            bits=np.float32(x).view(np.uint32).item()
            if bits == 0x7f800002:break
            result.append(None if bits == 0x7f800001 else float(x))
        else:
            raise ValidationError('Unsupported FORMAT type in benchmark adapter')
    return text(result[0] if scalar and result else (None if scalar else tuple(result)))


def carrier_calls(wrapper):
    record=wrapper.record
    # Unlike gt_types, the actual allele indexes retain partial ALT calls.
    gt=record.genotype.array()
    alleles=gt[:,:-1]
    # pysam 0.23.3 exposes out-of-range GT indexes as None on biallelic
    # records. Mirror that reader behavior, including partial 1/2 -> 1/. calls.
    alleles[alleles > 1] = -1
    if np.any(alleles < -2):
        raise ValidationError('Invalid genotype sentinel')
    selected=np.flatnonzero(np.any(alleles == 1,axis=1))
    fields={k:record.format(k) if k in wrapper.format else None for k in ['GQ','DP','AD','FT']}
    site_filter=';'.join(record.FILTERS) or '.'
    # cyvcf2 collapses phasing flags for polyploids differently from pysam.
    # Only this uncommon case needs textual GT separators to preserve the
    # production rule (all separators phased); diploid arrays stay vectorized.
    raw_samples=str(record).rstrip('\n').split('\t')[9:] if alleles.shape[1]>2 else None
    gt_column=wrapper.format.index('GT')
    for i in selected:
        call=[None if int(a)==-1 else int(a) for a in alleles[i] if int(a)!=-2]
        phased=bool(gt[i,-1])
        if raw_samples is not None:
            raw_gt=raw_samples[i].split(':')[gt_column]
            phased='|' in raw_gt and '/' not in raw_gt
        row=dict(CHROM=wrapper.contig,POS=wrapper.pos,REF=wrapper.ref,ALT=wrapper.alts[0],
                 sample=wrapper.reader.header.samples[i],
                 GT=('|' if phased else '/').join(text(a) for a in call),
                 alt_dosage=call.count(1),site_FILTER=site_filter)
        row.update({k:format_value(fields[k],i,wrapper.reader.types.get(k,{})) for k in fields})
        yield row,None in call


def install(production):
    """Patch only this benchmark child process; share production logic and BGZF writer."""
    original=production.pysam
    production.pysam=SimpleNamespace(VariantFile=Reader,BGZFile=original.BGZFile,__version__=original.__version__)
    production.carrier_calls=carrier_calls
