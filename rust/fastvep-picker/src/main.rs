use std::cmp::Ordering;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const OUTPUT_FIELDS: [&str; 21] = [
    "record_id",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "Feature",
    "Gene",
    "SYMBOL",
    "BIOTYPE",
    "Consequence",
    "EXON",
    "INTRON",
    "CDS_position",
    "Protein_position",
    "HGVS_OFFSET",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "CANONICAL",
    "TSL",
    "APPRIS",
    "CCDS",
];

#[derive(Debug)]
struct Args {
    fastvep: String,
    transcript_priority: PathBuf,
    consequence_ranks: PathBuf,
    output: PathBuf,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct Priority {
    mane_select: i64,
    mane_plus_clinical: i64,
    canonical: i64,
    appris: i64,
    tsl: i64,
    protein_coding: i64,
    ccds: i64,
    length: i64,
    ensembl: i64,
    refseq: i64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct PickKey {
    no_mane_select: i64,
    no_mane_plus_clinical: i64,
    no_canonical: i64,
    appris: i64,
    tsl: i64,
    non_protein_coding: i64,
    no_ccds: i64,
    consequence_rank: i64,
    length: i64,
    ensembl: i64,
    refseq: i64,
}

impl Ord for PickKey {
    fn cmp(&self, other: &Self) -> Ordering {
        (
            self.no_mane_select,
            self.no_mane_plus_clinical,
            self.no_canonical,
            self.appris,
            self.tsl,
            self.non_protein_coding,
            self.no_ccds,
            self.consequence_rank,
            self.length,
            self.ensembl,
            self.refseq,
        )
            .cmp(&(
                other.no_mane_select,
                other.no_mane_plus_clinical,
                other.no_canonical,
                other.appris,
                other.tsl,
                other.non_protein_coding,
                other.no_ccds,
                other.consequence_rank,
                other.length,
                other.ensembl,
                other.refseq,
            ))
    }
}

impl PartialOrd for PickKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn usage() -> &'static str {
    "Usage: fastvep-picker --fastvep <VCF|-> --transcript-priority <TSV> \\\n+     --consequence-ranks <TSV> --output <TSV>"
}

fn parse_args() -> Result<Args, String> {
    let mut values: HashMap<String, String> = HashMap::new();
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            println!("{}", usage());
            std::process::exit(0);
        }
        if !arg.starts_with("--") {
            return Err(format!("unexpected argument: {arg}\n{}", usage()));
        }
        let value = args
            .next()
            .ok_or_else(|| format!("missing value for {arg}\n{}", usage()))?;
        values.insert(arg, value);
    }
    let take = |name: &str| {
        values
            .get(name)
            .cloned()
            .ok_or_else(|| format!("missing {name}\n{}", usage()))
    };
    Ok(Args {
        fastvep: take("--fastvep")?,
        transcript_priority: PathBuf::from(take("--transcript-priority")?),
        consequence_ranks: PathBuf::from(take("--consequence-ranks")?),
        output: PathBuf::from(take("--output")?),
    })
}

fn parse_i64(value: Option<&&str>, default: i64) -> i64 {
    value.and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn read_tsv(path: &Path) -> io::Result<(Vec<String>, Vec<Vec<String>>)> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty TSV"))?
        .split('\t')
        .map(str::to_owned)
        .collect();
    let rows = lines
        .map(|line| line.map(|text| text.split('\t').map(str::to_owned).collect::<Vec<_>>()))
        .collect::<io::Result<Vec<_>>>()?;
    Ok((header, rows))
}

fn field_index(header: &[String], name: &str) -> io::Result<usize> {
    header
        .iter()
        .position(|field| field == name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing TSV field {name}"),
            )
        })
}

fn load_priorities(path: &Path) -> io::Result<HashMap<String, Priority>> {
    let (header, rows) = read_tsv(path)?;
    let indexes = [
        "transcript_id",
        "mane_select",
        "mane_plus_clinical",
        "canonical",
        "appris",
        "tsl",
        "protein_coding",
        "ccds",
        "length",
        "ensembl",
        "refseq",
    ]
    .map(|name| field_index(&header, name))
    .into_iter()
    .collect::<io::Result<Vec<_>>>()?;
    let mut result = HashMap::with_capacity(rows.len());
    for row in rows {
        let values: Vec<&str> = row.iter().map(String::as_str).collect();
        let at = |index: usize| values.get(index);
        result.insert(
            at(indexes[0]).copied().unwrap_or_default().to_owned(),
            Priority {
                mane_select: parse_i64(at(indexes[1]), 0),
                mane_plus_clinical: parse_i64(at(indexes[2]), 0),
                canonical: parse_i64(at(indexes[3]), 0),
                appris: parse_i64(at(indexes[4]), 100),
                tsl: parse_i64(at(indexes[5]), 100),
                protein_coding: parse_i64(at(indexes[6]), 0),
                ccds: parse_i64(at(indexes[7]), 0),
                length: parse_i64(at(indexes[8]), 0),
                ensembl: parse_i64(at(indexes[9]), 1),
                refseq: parse_i64(at(indexes[10]), 1),
            },
        );
    }
    Ok(result)
}

fn load_ranks(path: &Path) -> io::Result<HashMap<String, i64>> {
    let (header, rows) = read_tsv(path)?;
    let consequence = field_index(&header, "consequence")?;
    let rank = field_index(&header, "rank")?;
    Ok(rows
        .into_iter()
        .filter_map(|row| Some((row.get(consequence)?.clone(), row.get(rank)?.parse().ok()?)))
        .collect())
}

fn bare_transcript(value: &str) -> &str {
    value.split_once('.').map_or(value, |(bare, _)| bare)
}

fn enabled(value: i64) -> i64 {
    if value != 0 { 0 } else { 1 }
}

fn build_key(
    values: &[&str],
    feature_index: usize,
    consequence_index: usize,
    priorities: &HashMap<String, Priority>,
    ranks: &HashMap<String, i64>,
) -> PickKey {
    let feature = values.get(feature_index).copied().unwrap_or_default();
    let default_priority = Priority {
        appris: 100,
        tsl: 100,
        ensembl: 1,
        refseq: 1,
        ..Priority::default()
    };
    let priority = priorities
        .get(bare_transcript(feature))
        .unwrap_or(&default_priority);
    let consequence_rank = values
        .get(consequence_index)
        .copied()
        .unwrap_or_default()
        .split('&')
        .map(|term| ranks.get(term).copied().unwrap_or(1000))
        .min()
        .unwrap_or(1000);
    PickKey {
        no_mane_select: enabled(priority.mane_select),
        no_mane_plus_clinical: enabled(priority.mane_plus_clinical),
        no_canonical: enabled(priority.canonical),
        appris: priority.appris,
        tsl: priority.tsl,
        non_protein_coding: enabled(priority.protein_coding),
        no_ccds: enabled(priority.ccds),
        consequence_rank,
        length: priority.length,
        ensembl: priority.ensembl,
        refseq: priority.refseq,
    }
}

fn input_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    if path == "-" {
        Ok(Box::new(BufReader::with_capacity(1024 * 1024, io::stdin())))
    } else {
        Ok(Box::new(BufReader::with_capacity(
            1024 * 1024,
            File::open(path)?,
        )))
    }
}

fn csq_schema(line: &str) -> Option<Vec<String>> {
    if !line.starts_with("##INFO=<ID=CSQ") {
        return None;
    }
    let start = line.find("Format: ")? + "Format: ".len();
    let rest = &line[start..];
    let end = rest.find('"').unwrap_or(rest.len());
    Some(rest[..end].split('|').map(str::to_owned).collect())
}

fn process(args: &Args) -> io::Result<u64> {
    let priorities = load_priorities(&args.transcript_priority)?;
    let ranks = load_ranks(&args.consequence_ranks)?;
    let mut reader = input_reader(&args.fastvep)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(&args.output)?);
    writeln!(writer, "{}", OUTPUT_FIELDS.join("\t"))?;

    let mut schema: Option<Vec<String>> = None;
    let mut schema_indexes = Vec::new();
    let mut feature_index = 0;
    let mut consequence_index = 0;
    let mut line = String::with_capacity(4096);
    let mut records = 0_u64;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with("##INFO=<ID=CSQ") {
            let fields = csq_schema(&line).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "could not parse CSQ schema")
            })?;
            schema_indexes = OUTPUT_FIELDS[5..]
                .iter()
                .map(|name| {
                    fields
                        .iter()
                        .position(|field| field == name)
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("CSQ schema missing field {name}"),
                            )
                        })
                })
                .collect::<io::Result<Vec<_>>>()?;
            feature_index = fields.iter().position(|field| field == "Feature").unwrap();
            consequence_index = fields
                .iter()
                .position(|field| field == "Consequence")
                .unwrap();
            schema = Some(fields);
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        if schema.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no CSQ schema found before records",
            ));
        }
        let columns: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
        if columns.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "VCF row has fewer than 8 columns",
            ));
        }
        let csq = columns[7]
            .split(';')
            .find_map(|value| value.strip_prefix("CSQ="));
        let Some(csq) = csq else { continue };

        let mut selected: Option<Vec<&str>> = None;
        let mut selected_key: Option<PickKey> = None;
        for annotation in csq.split(',') {
            let values: Vec<&str> = annotation.split('|').collect();
            let key = build_key(
                &values,
                feature_index,
                consequence_index,
                &priorities,
                &ranks,
            );
            if selected_key.as_ref().is_none_or(|current| key < *current) {
                selected_key = Some(key);
                selected = Some(values);
            }
        }
        let Some(selected) = selected else { continue };
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            columns[2], columns[0], columns[1], columns[3], columns[4]
        )?;
        for index in &schema_indexes {
            write!(
                writer,
                "\t{}",
                selected.get(*index).copied().unwrap_or_default()
            )?;
        }
        writeln!(writer)?;
        records += 1;
    }
    writer.flush()?;
    Ok(records)
}

fn main() {
    let args = parse_args().unwrap_or_else(|message| {
        eprintln!("{message}");
        std::process::exit(2);
    });
    match process(&args) {
        Ok(records) => println!("records={records}"),
        Err(error) => {
            eprintln!("fastvep-picker: {error}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bare_transcript_removes_version() {
        assert_eq!(bare_transcript("ENST0001.4"), "ENST0001");
        assert_eq!(bare_transcript("ENST0001"), "ENST0001");
    }

    #[test]
    fn mane_beats_more_severe_consequence() {
        let mut priorities = HashMap::new();
        priorities.insert(
            "MANE".to_owned(),
            Priority {
                mane_select: 1,
                ..Priority::default()
            },
        );
        let ranks = HashMap::from([
            ("stop_gained".to_owned(), 1),
            ("missense_variant".to_owned(), 10),
        ]);
        let other = vec!["stop_gained", "OTHER"];
        let mane = vec!["missense_variant", "MANE"];
        assert!(
            build_key(&mane, 1, 0, &priorities, &ranks)
                < build_key(&other, 1, 0, &priorities, &ranks)
        );
    }
}
