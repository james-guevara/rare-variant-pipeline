# GitHub Preparation Checklist

## Completed ✓

- [x] Fixed `family_query.py` script (removed debug code, completed implementation)
- [x] Added comprehensive `README.md` documentation
- [x] Created `.gitignore` file
- [x] Added docstrings to all Python scripts
- [x] Reviewed all scripts for code quality

## Before First Commit

### 1. Clean Up Large/Temporary Files

**Core dumps (~27GB total) - RECOMMENDED TO DELETE:**
```bash
rm core.2246369 core.2246748 core.2247352 core.2247533 core.2248210 core.2248578
```

**Test/temporary files:**
```bash
rm tmp.tsv duplicated_rows.tsv test.variants.tsv
```

**Optional - Clean up logs (can keep for debugging):**
```bash
# Only if you don't need these logs
# rm -rf logs_*/*.out logs_*/*.err
```

### 2. Review Symlinks

These symlinks point to external resources and won't transfer to GitHub:
- `out_vep` → external path
- `resources` → external path
- `vcfs` → external path
- `VEP_CACHE` → external path
- `VEP_PLUGINS` → external path
- `VEP_PLUGINS_ALL` → external path

**Action:** Document the expected structure in README (already done)

### 3. Verify What Will Be Committed

```bash
# Check current status
git status

# See what will be tracked
git ls-files

# Check file sizes
du -sh *
```

### 4. Update README with Your Information

Edit `README.md` and fill in:
- **License** section (line ~270): Add your preferred license (MIT, GPL-3.0, etc.)
- **Contact** section (line ~272): Add your contact information
- Any project-specific details (institution, funding, etc.)

## Initialize Git Repository

```bash
# Initialize repository
git init

# Add files
git add .gitignore
git add README.md
git add *.py
git add *.sh
git add GITHUB_PREP.md  # optional

# First commit
git commit -m "Initial commit: Rare variant annotation pipeline

- VEP-based variant annotation workflow
- Family-based genotype extraction
- Multi-allelic variant resolution
- Constraint metric integration
- Complete documentation"
```

## Create GitHub Repository

### Option 1: Using GitHub CLI (gh)

```bash
# Create private repository
gh repo create rare-variant-pipeline --private --source=. --remote=origin

# Or create public repository
gh repo create rare-variant-pipeline --public --source=. --remote=origin

# Push to GitHub
git push -u origin main
```

### Option 2: Using Web Interface

1. Go to https://github.com/new
2. Create repository named `rare-variant-pipeline` (or your preferred name)
3. Choose public or private
4. **Do NOT** initialize with README, .gitignore, or license
5. Follow the instructions to push existing repository:

```bash
git remote add origin https://github.com/YOUR_USERNAME/rare-variant-pipeline.git
git branch -M main
git push -u origin main
```

## Recommended Additional Steps

### 1. Add a LICENSE File

Example for MIT License:
```bash
cat > LICENSE << 'EOF'
MIT License

Copyright (c) 2024 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF
```

### 2. Add Requirements File

```bash
cat > requirements.txt << 'EOF'
polars>=0.19.0
polars-bio>=0.1.0
cyvcf2>=0.30.0
numpy>=1.24.0
EOF
```

### 3. Add GitHub Actions (Optional)

For automated testing, create `.github/workflows/test.yml`

### 4. Add Contributing Guidelines (Optional)

Create `CONTRIBUTING.md` if you want to accept contributions

## Repository Structure After Setup

```
rare-variant-pipeline/
├── .git/
├── .gitignore
├── README.md
├── LICENSE
├── requirements.txt
├── GITHUB_PREP.md
├── RUN_ANNOTATION_PIPELINE.sb.sh
├── RUN_ANNOTATION_PIPELINE.postprocess.sh
├── RUN_ANNOTATION_PIPELINE.family_query.sb.sh
├── RUN_ANNOTATION_PIPELINE.genotype_query.sb.sh
├── RUN_ANNOTATION_PIPELINE.make_tables.sh
├── RUN_ANNOTATION_PIPELINE.make_final_table.sh
├── family_query.py
├── reformat_variants.py
├── resolve_genotypes.py
├── resolve_family_genotypes.py
├── merge_family_variants.py
├── concat_family_merged.py
├── containers/              # Container definitions
├── tests/                   # Test scripts
└── wgs.psam                 # Sample metadata
```

## Files Excluded by .gitignore

The following will NOT be committed (good!):
- Core dumps (`core.*`)
- Output directories (`out/`, `out_vep/`, `logs_*/`)
- Large data files (`*.vcf.gz`, `*.tsv.gz`)
- Python cache (`__pycache__/`, `*.pyc`)
- Temporary files (`tmp/`, `*.tmp`)
- SLURM logs (`*.out`, `*.err`)

## Troubleshooting

### If repository is too large:
```bash
# Check what's taking up space
du -sh * | sort -h

# Remove from git if accidentally added
git rm --cached <large-file>
git commit --amend
```

### If you need to add more to .gitignore after first commit:
```bash
# Edit .gitignore
nano .gitignore

# Remove cached files
git rm -r --cached .

# Re-add everything (respecting new .gitignore)
git add .
git commit -m "Update .gitignore"
```

## Next Steps After GitHub Setup

1. Add repository badge to README
2. Create releases/tags for versions
3. Set up GitHub Issues for bug tracking
4. Add GitHub wiki for extended documentation
5. Consider adding example datasets (small test data)
6. Add Docker/Singularity container recipes

## Questions?

- GitHub Docs: https://docs.github.com/
- Git Basics: https://git-scm.com/book/en/v2
- Licensing Help: https://choosealicense.com/
