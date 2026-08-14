#!/usr/bin/env python3
"""
extract_orthogroup_gene_info.py

For a given list of orthogroups, look up the protein(s) for each
species, then trace that protein through the species' GFF3 annotation back to
its parent gene, and report the gene's genomic coordinates and name. Output is a 
TSV with one row per (orthogroup, species, protein):

USAGE: 
#   python3 extract_orthogroup_gene_info.py Orthogroups_SingleCopyOrthologues.txt \ # input file name can be changed

EXPLICIT USAGE:

python3 extract_orthogroup_gene_info.py \
    --orthogroups-tsv Orthogroups.tsv \
    --gff-dir 'gff3_files' \ # directory containing the <species>.annotation.gff3 files
    --output orthogroup_gene_info.tsv # output file name




The lookup chain for one protein is:

    protein ID  --(CDS row in GFF3)-->  GeneID  --(gene row in GFF3)-->  gene info

    CDS row:   ...;Dbxref=GeneID:101805557,GenBank:...;Name=NP_001265621.1;...
    gene row:  ...;Dbxref=GeneID:101805557;Name=rbpja;...   (cols give chrom/start/end/strand)

Output is a TSV with one row per (orthogroup, species, protein):
    orthogroup, species, protein_id, gene_id, gene_name, chromosome, start, end, strand


FUNCTIONS
# read_orthogroup_ids           - reads the input file and returns a list of orthogroup IDs
# parse_orthogroups_tsv         - reads Orthogroups.tsv and maps each orthogroup to its protein IDs per species
# parse_gff_attributes          - converts a GFF3 attributes string (col 9) into a dictionary
# extract_gene_id               - pulls the numeric GeneID out of the Dbxref attribute
# scan_gff3                     - reads a GFF3 file once and returns two lookups: protein->GeneID and GeneID->coordinates
# collect_proteins_per_species  - flattens ortho_to_proteins into a per-species set of all proteins needed
# resolve_gene_info             - for each species, scans its GFF3 and resolves every protein to its gene coordinates
# write_output                  - writes the final TSV file with one row per orthogroup/species/protein
# main                          - reads command-line arguments and runs all the above functions in order


Design notes
------------
The GFF3 files are very large, so scan each file
exactly once, line by line, instead of re-searching the file for every protein.
During that single pass we collect:
  * for the proteins we actually care about: protein_id -> GeneID  (from CDS rows)
  * for every gene: GeneID -> (chrom, start, end, strand, gene_name)  (from gene rows)
and then join the two afterwards.
"""

import argparse
import os
import sys


# ---------------------------------------------------------------------------
# Step 1: Read the list of orthogroup IDs we want to process.
# ---------------------------------------------------------------------------
def read_orthogroup_ids(path):
    """Read a plain-text file of orthogroup IDs (one per line, no header).

    Blank lines and surrounding whitespace are ignored. Returns the IDs as a
    list, preserving the order in which they appear in the file so the output
    is reproducible.
    """
    orthogroup_ids = []
    with open(path) as handle:
        for line in handle:
            ortho_id = line.strip()
            if ortho_id:  # skip empty / whitespace-only lines
                orthogroup_ids.append(ortho_id)
    return orthogroup_ids


# ---------------------------------------------------------------------------
# Step 2: Parse Orthogroups.tsv into a lookup of orthogroup -> species -> proteins.
# ---------------------------------------------------------------------------
def parse_orthogroups_tsv(path, wanted_ids):
    """Parse Orthogroups.tsv, keeping only the orthogroups we asked for.

    The TSV has a header whose first column is "Orthogroup" and whose remaining
    columns are named like "Metriaclima_zebra_AF.protein". Each data cell holds
    a comma-separated list of protein IDs for that species (for single-copy
    orthogroups this is normally just one protein).

    Returns a tuple of:
      * species_names: list of species base names (header minus ".protein"),
        in column order.
      * ortho_to_proteins: dict mapping orthogroup ID ->
        { species_name: [protein_id, ...] }.

    `wanted_ids` is used as a filter so we don't carry orthogroups we don't need.
    """
    wanted = set(wanted_ids)
    ortho_to_proteins = {}

    with open(path) as handle:
        # The header tells us which column belongs to which species.
        header = handle.readline().rstrip("\n").split("\t")
        # Drop the ".protein" suffix so column "Metriaclima_zebra_AF.protein"
        # becomes the species base name "Metriaclima_zebra_AF". Skip column 0
        # ("Orthogroup"), which is not a species.
        species_names = [col[:-len(".protein")] if col.endswith(".protein") else col
                         for col in header[1:]]

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            ortho_id = fields[0]
            if ortho_id not in wanted:
                continue

            # Pair each species with its cell value. Cells may be missing at the
            # end of a line, so we guard against a short row.
            species_to_proteins = {}
            for index, species in enumerate(species_names):
                cell = fields[index + 1] if (index + 1) < len(fields) else ""
                # Split the comma-separated protein list and trim whitespace.
                proteins = [p.strip() for p in cell.split(",") if p.strip()]
                species_to_proteins[species] = proteins

            ortho_to_proteins[ortho_id] = species_to_proteins

    return species_names, ortho_to_proteins


# ---------------------------------------------------------------------------
# Step 3: Helpers for reading the GFF3 attributes column (column 9).
# ---------------------------------------------------------------------------
def parse_gff_attributes(attribute_field):
    """Turn a GFF3 attributes string into a dict.

    The attributes column looks like:
        ID=cds-XP_1;Parent=rna-XM_1;Dbxref=GeneID:123,GenBank:XP_1;Name=foo
    i.e. semicolon-separated key=value pairs. Returns {"ID": "cds-XP_1", ...}.
    """
    attributes = {}
    for pair in attribute_field.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            attributes[key] = value
    return attributes


def extract_gene_id(attributes):
    """Pull the numeric GeneID out of a parsed attributes dict.

    The Dbxref value bundles several cross-references separated by commas, e.g.
        GeneID:101805557,GenBank:NP_001265621.1
    We want only the "GeneID:" token, and we return just the number after it.
    Returns None if no GeneID is present.
    """
    dbxref = attributes.get("Dbxref", "")
    for token in dbxref.split(","):
        if token.startswith("GeneID:"):
            return token[len("GeneID:"):]
    return None


# ---------------------------------------------------------------------------
# Step 4: Single-pass scan of one species' GFF3 file.
# ---------------------------------------------------------------------------
def scan_gff3(path, wanted_protein_ids):
    """Scan one GFF3 file once and gather everything we need from it.

    Two things are collected in the same pass:
      * protein_to_geneid: for CDS rows whose Name is a protein we care about,
        map protein_id -> GeneID.
      * geneid_to_info: for every gene row, map GeneID ->
        (chromosome, start, end, strand, gene_name).

    Both maps are returned so the caller can join protein -> GeneID -> gene info.
    A gene row may appear before or after its CDS rows, which is why we simply
    record all genes rather than trying to look them up on the fly.
    """
    wanted_protein_ids = set(wanted_protein_ids)
    protein_to_geneid = {}
    geneid_to_info = {}

    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):  # comment / header lines
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:  # not a valid feature line
                continue

            feature_type = fields[2]

            # We only care about two feature types: CDS and gene.
            if feature_type == "CDS":
                attributes = parse_gff_attributes(fields[8])
                # The CDS's Name is the protein ID; only keep ones we need.
                protein_id = attributes.get("Name")
                if protein_id in wanted_protein_ids:
                    gene_id = extract_gene_id(attributes)
                    if gene_id is not None:
                        # A protein has many CDS rows; they all share the same
                        # GeneID, so recording it once (or repeatedly) is fine.
                        protein_to_geneid[protein_id] = gene_id

            elif feature_type == "gene":
                attributes = parse_gff_attributes(fields[8])
                gene_id = extract_gene_id(attributes)
                if gene_id is not None:
                    geneid_to_info[gene_id] = (
                        fields[0],               # chromosome / sequence ID
                        fields[3],               # start
                        fields[4],               # end
                        fields[6],               # strand
                        attributes.get("Name", ""),  # human-readable gene name
                    )

    return protein_to_geneid, geneid_to_info


# ---------------------------------------------------------------------------
# Step 5: Figure out which proteins each species' GFF3 needs to find.
# ---------------------------------------------------------------------------
def collect_proteins_per_species(ortho_to_proteins):
    """Build species_name -> set of all protein IDs needed for that species.

    Gathering the full set up front lets us scan each GFF3 only once and filter
    for exactly the proteins of interest.
    """
    proteins_per_species = {}
    for species_to_proteins in ortho_to_proteins.values():
        for species, proteins in species_to_proteins.items():
            proteins_per_species.setdefault(species, set()).update(proteins)
    return proteins_per_species


# ---------------------------------------------------------------------------
# Step 6: Resolve every protein to its gene info, species by species.
# ---------------------------------------------------------------------------
def resolve_gene_info(species_names, proteins_per_species, gff_dir):
    """For each species, scan its GFF3 and resolve protein -> gene info.

    Returns a nested dict: species_name -> protein_id -> info-dict, where each
    info-dict holds gene_id, gene_name, chromosome, start, end and strand.
    Proteins that can't be resolved are simply absent (the caller decides how to
    report the gap).
    """
    resolved = {}

    for species in species_names:
        wanted_proteins = proteins_per_species.get(species, set())
        if not wanted_proteins:
            resolved[species] = {}
            continue

        # GFF3 filenames follow the pattern "<species>.annotation.gff3".
        gff_path = os.path.join(gff_dir, species + ".annotation.gff3")
        if not os.path.isfile(gff_path):
            # Warn but keep going so the other species still get processed.
            print("WARNING: GFF3 not found for species '%s': %s"
                  % (species, gff_path), file=sys.stderr)
            resolved[species] = {}
            continue

        print("Scanning %s (%d proteins to find)..."
              % (os.path.basename(gff_path), len(wanted_proteins)),
              file=sys.stderr)

        protein_to_geneid, geneid_to_info = scan_gff3(gff_path, wanted_proteins)

        # Join the two maps: protein -> GeneID -> gene info.
        species_resolved = {}
        for protein_id, gene_id in protein_to_geneid.items():
            info = geneid_to_info.get(gene_id)
            if info is None:
                # We found the CDS but not a matching gene row.
                continue
            chromosome, start, end, strand, gene_name = info
            species_resolved[protein_id] = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "strand": strand,
            }
        resolved[species] = species_resolved

    return resolved


# ---------------------------------------------------------------------------
# Step 7: Write the final TSV, assembling one row per (orthogroup, species, protein).
# ---------------------------------------------------------------------------
def write_output(output_path, orthogroup_ids, species_names,
                 ortho_to_proteins, resolved):
    """Write the results to a TSV file.

    We iterate orthogroups in their original input order, and species in their
    TSV column order, so the output is deterministic. Any protein that could
    not be resolved to a gene is still emitted, with empty gene fields, so the
    user can see which lookups failed rather than silently dropping them.
    """
    columns = ["orthogroup", "species", "protein_id", "gene_id",
               "gene_name", "chromosome", "start", "end", "strand"]

    rows_written = 0
    with open(output_path, "w") as out:
        out.write("\t".join(columns) + "\n")

        for ortho_id in orthogroup_ids:
            species_to_proteins = ortho_to_proteins.get(ortho_id)
            if species_to_proteins is None:
                # The orthogroup ID wasn't found in Orthogroups.tsv at all.
                print("WARNING: orthogroup '%s' not found in Orthogroups.tsv"
                      % ortho_id, file=sys.stderr)
                continue

            for species in species_names:
                for protein_id in species_to_proteins.get(species, []):
                    info = resolved.get(species, {}).get(protein_id)
                    if info is None:
                        # Protein could not be traced to a gene; emit blanks.
                        info = {"gene_id": "", "gene_name": "", "chromosome": "",
                                "start": "", "end": "", "strand": ""}

                    out.write("\t".join([
                        ortho_id,
                        species,
                        protein_id,
                        info["gene_id"],
                        info["gene_name"],
                        info["chromosome"],
                        info["start"],
                        info["end"],
                        info["strand"],
                    ]) + "\n")
                    rows_written += 1

    print("Wrote %d rows to %s" % (rows_written, output_path), file=sys.stderr)


# ---------------------------------------------------------------------------
# Step 8: Command-line interface and orchestration.
# ---------------------------------------------------------------------------
def main():
    """Parse arguments and run the full pipeline end to end."""
    parser = argparse.ArgumentParser(
        description="Map orthogroup single-copy proteins to their gene "
                    "coordinates using species GFF3 annotations.")
    parser.add_argument(
        "orthogroups_file",
        help="Input file of orthogroup IDs, one per line, no header "
             "(e.g. Orthogroups_SingleCopyOrthologues.txt).")
    parser.add_argument(
        "--orthogroups-tsv", default="Orthogroups.tsv",
        help="Orthogroups.tsv mapping orthogroups to per-species proteins "
             "(default: Orthogroups.tsv).")
    parser.add_argument(
        "--gff-dir", default=".",
        help="Directory containing the <species>.annotation.gff3 files "
             "(default: current directory).")
    parser.add_argument(
        "--output", default="orthogroup_gene_info.tsv",
        help="Path to write the output TSV "
             "(default: orthogroup_gene_info.tsv).")
    args = parser.parse_args()

    # Run the steps in order, passing results from one stage to the next.
    orthogroup_ids = read_orthogroup_ids(args.orthogroups_file)
    species_names, ortho_to_proteins = parse_orthogroups_tsv(
        args.orthogroups_tsv, orthogroup_ids)
    proteins_per_species = collect_proteins_per_species(ortho_to_proteins)
    resolved = resolve_gene_info(species_names, proteins_per_species, args.gff_dir)
    write_output(args.output, orthogroup_ids, species_names,
                 ortho_to_proteins, resolved)


if __name__ == "__main__":
    main()
