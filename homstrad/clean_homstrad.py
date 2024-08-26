"""
Splits Homstard superposed PDB files (*-sup.pdb) into individual PDB files for structure aligners.
Also generates an amino acid FASTA file (for sequence aligners).

python3 clean_homstrad.py homstrad/<family> outputfolder/
"""

import argparse
import re

from pathlib import Path
from collections import defaultdict


def main(folder_path, output_path):
    # Inputs
    folder = Path(folder_path)
    alignment = folder / f"{folder.stem}.ali"
    superposition = folder / f"{folder.stem}-sup.pdb"
    assert folder.is_dir(), f"{folder} is not a directory"
    assert alignment.exists(), f"Could not find {alignment}"
    assert superposition.exists(), f"Could not find {superposition}"
    print("Inputs", folder, alignment, superposition)

    output = Path(output_path)
    if not output.is_dir():
        print(f"Making directory {output}")
        output.mkdir(parents=True)

    output_pdbs = output / "pdbs"
    if not output_pdbs.is_dir():
        print(f"Making directory {output_pdbs}")
        output_pdbs.mkdir(parents=True)

    # Get family members directly from .ali file
    # Matches blocks from > up until sequence terminating *
    with alignment.open() as fp:
        text = fp.read()
    members = []
    sequences = []
    for name, sequence in re.findall(r"^>.+?;(?P<name>.+?)$\nstructure[A-Z]?.+?$\n(?P<sequence>.+?)\*", text, re.MULTILINE | re.DOTALL):
        members.append(name)
        sequences.append(str(sequence).replace('\n', '').replace('/', '-'))

    # Read superposition PDB file
    with superposition.open() as fp:
        text = fp.read()

    # Find chain assignments in REMARK section, if present e.g.
    # REMARK   1bif   chain A
    # REMARK   1k6ma  chain B
    # --> { 'A': '1bif', 'B': '1k6ma' }
    remark = {
        chain: name
        for (name, chain) in
        re.findall(r"^REMARK\s*(?P<name>\w+?)\s*chain\s*(?P<chain>[A-Z0-9]*)\s*?$", text, re.MULTILINE)
    }
    
    # Split the superposed PDB into individual files
    # ATOM lines saved into explicitly mapped chain if remarks found,
    # otherwise just following .ali file order
    single_pdbs = defaultdict(str)
    last_chain = None
    last_residue = None
    chain_index = 0
    residue_index = 0
    name = ""
    for line in text.split('\n'):
        if line.startswith("ATOM"):
            chain = line[21:22]
            residue = line[22:27]
            if residue != last_residue:
                residue_index += 1
                last_residue = residue
            if chain != last_chain:
                name = remark.get(chain, members[chain_index])
                chain_index += 1
                residue_index = 0
                last_chain = chain
            line = line[0:22] + f"{residue_index:-4}" + line[26:] + '\n'
            single_pdbs[name] += line

    # Write out member.pdb files into output folder
    for member in members:
        assert member in single_pdbs, f"{member} not parsed"
        pdb = output_pdbs / f"{member}.pdb"
        print(f"Writing {pdb}")
        with pdb.open('w') as fp:
            fp.write(single_pdbs[member])

    # Write Homstrad MSA as clean FASTA file
    msa_output = output / f"{folder.stem}_msa.fasta"
    print(f"Writing {msa_output}")
    with msa_output.open('w') as fp:
        records = "".join(f">{name}\n{sequence}\n" for name, sequence in zip(members, sequences))
        fp.write(records)

    # Filter out gaps for raw AA sequence
    aa_output = output / f"{folder.stem}_aa.fasta"
    print(f"Writing {aa_output}")
    with aa_output.open('w') as fp:
        records = "".join(f">{name}\n{sequence.replace('-', '')}\n" for name, sequence in zip(members, sequences))
        fp.write(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split Homstrad *-sup.pdb into distinct PDB files and write AA FASTA file")
    parser.add_argument("folder_path", help="Homstrad family folder")
    parser.add_argument("output_path", help="Output folder")
    args = parser.parse_args()
    main(args.folder_path, args.output_path)
