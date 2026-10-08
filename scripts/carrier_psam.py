"""Sample metadata only: no sex interpretation, frequency counting or subsetting."""
import csv
from extract_exact_carriers import ValidationError

FIELDS = ['#IID', 'SEX', 'participant_id', 'frequency_representative', 'unrelated']


def load_psam(path, samples):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(FIELDS) <= set(header):
            raise ValidationError('PSAM requires unique columns: #IID SEX participant_id frequency_representative unrelated')
        rows = {}
        for row in reader:
            if None in row or any(row.get(k) is None for k in header):
                raise ValidationError('Malformed PSAM row')
            if any(not row[k].strip() for k in FIELDS):
                raise ValidationError('Empty required PSAM value; use explicit unknown SEX encoding')
            if row['#IID'] in rows:
                raise ValidationError('Duplicate PSAM sample identifier')
            for field in ['frequency_representative', 'unrelated']:
                if row[field] not in ('0', '1'):
                    raise ValidationError('PSAM frequency flags must be explicit 0 or 1')
            rows[row['#IID']] = row
    if not set(samples) <= rows.keys():
        raise ValidationError('PSAM missing source VCF samples')
    # Keep VCF order, not PSAM order. Extra metadata samples need not be in a block.
    return header, [rows[s] for s in samples], len(rows) - len(samples)
