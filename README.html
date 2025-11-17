Open Babel
obabel -ismi retinol.smi -oPDB --gen3d -O retinol.pdb
obabel -ismi retinol.smi -oPDB --gen3d --partialcharge gasteiger -p 7.4 -O retinol.pdbqt
obabel -ipdb retinol.pdb -osmi -O retinol_pdb.smi
obabel -ipdbqt input.pdbqt -osmi -O output.smi
Cut -f1 output.smi > clean.smi

Comparing Protein Structures
Comparing 2 structures in pymol
Load 1n8k and 1ye3
File -> Get PDB -> 1n8k
File -> Get PDB -> 1ye3
Align 1ye3 to 1n8k using default parameters
1ye3 : Action -> Align -> To molecule -> 1n8k (Aligned to chain A)
Want to align 1ye3 to chain B of 1n8k
sele 1n8k_B, (chain B)
1ye3 : Action -> Align -> To selection -> 1n8k_B

Comparing 2 structures in chimera
Load 1n8k and 1ye3
File -> Fetch by ID -> 1n8k
File -> Fetch by ID -> 1ye3
Align 1ye3 to 1n8k using default parameters
Tools -> Structure Comparison -> Match maker
1n8k Reference
1ye3 Structure
Want to align 1ye3 to chain B of 1n8k
Tools -> Structure Comparison -> Match maker -> Specific chain in reference structure with best aligning chain in match structure -> 1n8k chain B

for LSQMAN pdb shall be clean: (do only if gives error)
adding missing chain identifier
obabel input.pdb -O fixed.pdb --addchain A
Fixing alignment
obabel input.pdb -O fixed.pdb –fix
Remove junk
grep '^ATOM\|^HETATM\|^END' input.pdb > cleaned.pdb
echo "END" >> cleaned.pdb

LSQMAN
re m1 1n8k.pdb
re m2 1ye3.pdb
Always renames the chain in alphabetical order
ex
Mol 1? m1
Range 1? Starting position-end position
Mol 2? m2
Range 2? Starting position
ex m1 A1-400 m2 A1 -> format
at ma
ex m1 A1-400 m2 A1
at all
ex m1 A1-400 m2 A1
Apply rotation matrix
ap m1 m2 -> m2 is moved upon m1
Write the object onto output file
wr m2 m2.pdb
chimera 1n8k.pdb m2.pdb
RMSD = Relative RMSD × Estimated RMSD for 2 random proteins

PDB Cleaning
Using Command-line:
Remove Anisou
grep -vi aniso 4PPY.pdb > noaniso.pdb
Remove Alternate conformers 
awk '((substr($0,17,1)=="A" || substr($0,17,1)==" ") && $1~"^ATOM$") {print substr($0,1,16), substr($0,18,63)}' noaniso.pdb > step2_noaltconf.pdb
Keep Only a particular chain
awk '($5~"^B$"){print $0}' step2_noaltconf.pdb > step3_chB.pdb
Remove Hydrogens
awk '($12!~"^H$"){print $0}' step3.pdb > step4.pdb
Remove HETATM
grep '^ATOM' step3_chainB.pdb > final_clean.pdb
To see chain identifiers
grep "^ATOM" yourfile.pdb | awk '{print substr($0,22,1)}' | sort | uniq

To see number of chains present
grep "^ATOM" yourfile.pdb | awk '{print substr($0,22,1)}' | sort | uniq | wc -l

Dealing with non standard residues:
open in chimera → select residue → look at non std as well as list of standard
once found for all nonstd residues do
awk '$0 ~ /^HETATM/ {print $4, $5}' 4gwb.pdb | grep -v -E "ALA|CYS|ASP|GLU|PHE|GLY|HIS|ILE|LYS|LEU|MET|ASN|PRO|GLN|ARG|SER|THR|VAL|TRP|TYR" | sort | uniq -c
residue number can be seen → modify in pymol
in modres record you can see what actual aminoacid we need to modify into
wizard → mutagenesis
file → export molecule → save → file format pdb → name
sed just edits name enough for lsqman and all but not for modelling, simulations and docking for that you need to explicitely correct
sed -i 's/ MSE / MET /g' file.pdb

Install:
conda install -c conda-forge pdbfixer
conda install -c bioconda pdb-tools
using pdb-tools:
Delete Heteroatoms
pdb_delhetatm 4GWB.pdb > step1_nohet.pdb
Delete Altconf
pdb_selaltloc -A step1_nohet.pdb > noaltcon.pdb
pdb_delelem -H 4GWB.pdb > noH.pdb
Residue re numbering
pdb_reres 4GWB.pdb > renum.odb
renaming to chain A
pdb_chain -A renum.pdb > chA.pdb
Delete anisou
pdb_delrecords ANISOU 4GWB.pdb > 4GWB_noanisou.pdb
pdb_delrecords ANISOU,HETATM,REMARK input.pdb > clean.pdb
Select a particular chain
pdb_selchain -A input.pdb > chainA.pdb
pdb_selchain -A,B input.pdb > chains_AB.pdb

using pdb-fixer:
Fix missing atoms:
pdbfixer input.pdb --add-atoms=heavy --output fixed_missing_atoms.pdb
Fix non standard residues: (automatically detects modres field and modifies)
pdbfixer input.pdb --replace-nonstandard --output stdres.pdb
Add missing loops:
pdbfixer input.pdb --add-residues --output addloops.pdb
Docking:
pdbfixer input.pdb \
--replace-nonstandard \
--remove-heterogens=all \
--add-atoms=heavy \
--output dock_ready.pdb
MD:
pdbfixer input.pdb \
--replace-nonstandard \
--add-residues \
--add-atoms=heavy \
--add-hydrogens --ph=7.0 \
--remove-heterogens=water \
--output md_ready.pdb
Modeller:
pdbfixer template.pdb \
--replace-nonstandard \
--add-atoms=heavy \
--remove-heterogens=all \
--output modeller_ready.pdb
Energy Minimization:
clean PDB
pdb to gromacs:
gmx pdb2gmx -f step2_noaltconf.pdb -water spc -ignh -o model.gro -p model.top

select 6

made a box of 1nm (atleast 0.8 was mentioned, to avoid clashes due to smaller dimentions made 1nm box)
gmx editconf -f model.gro -bt cubic -d 1.0 -o box.gro

add water:
gmx solvate -cp box.gro -cs spc216.gro -p model.top -o solvated.gro

created tpr file for adding ions and neutralizing:
gmx grompp -f em.mdp -p model.top -c solvated.gro -o ion.tpr

add ions 200mM concentration
gmx genion -s ion.tpr -neutral -conc 0.2 -o ion.gro -p model.top

select 13 (SOL)

made em-steep.mdp file and added required parameters
gmx grompp -c ion.gro -p model.top -f em-steep.mdp -o em-steep.tpr
For warning like The largest distance between excluded atoms is 0.560 nm between atom 597 and 607 either do -maxwarn or in em.mdp give rcoulomb = 1.0 and rvdw     = 1.0 and rerun

gmx mdrun -v -deffnm em-steep

if asked to cascade outputs:
gmx grompp -f em1.mdp -c ion.gro -p model.top -o em1.tpr
gmx mdrun -v -deffnm em1
make changes in mdp
gmx grompp -f em2.mdp -c em1.gro -p model.top -o em2.tpr
gmx mdrun -v -deffnm em2

made em-cg.mdp file and added required parameters (integrator = cg
nsteps = 2000
emtol = 0.1
emstep = 0.01)
gmx grompp -f em-cg.mdp -c em-steep.gro -p model.top -o em-cg.tpr
gmx mdrun -v -deffnm em-cg

gmx trjconv -f em-cg.gro -pbc mol -o em-cg-model.pdb -s em-cg.tpr
select 1 (protien)
gmx trjconv -s em.tpr -f em.gro -o em_noPBC.gro -pbc mol -center

gmx rms -s distorted.gro -f em_noPBC.gro -o rmsd_allll.xvg #use this command for rmsd calculation
First select group 0 for system->All atoms
then select group for c-alpha atoms
then select group for backbone atoms
It will generate 3 files with calculated rmsd values


gmx trjconv -s topol.tpr -f traj.xtc -o traj_nojump.xtc -pbc nojump -center

gmx trjconv -f em-cg.gro -pbc mol -o em-cg-mol.pdb -s em-cg.tpr
Select a group: 1

gmx trjconv -f em-steep.gro -pbc mol -o em-steep-mol.pdb -s em-steep.tpr

(Analyze the enrgies
em-steep.log
- Copy the start step energies into text file
- Copy the last step energies into text file
- energy has reduced a lot, but which component?
- LJ 14 -> long range Lennard Jones

em-cg.log
- Copy the start step energies into text file
- Copy the last step energies into text file
- final energy has not changed much)

gmx energy -f em-steep.edr

xmgrace energy.xvg
(data -> import -> Filter: replace /*.dat with /energy.xvg -> Enter -> Load as NXY -> Apply)
Overlay plots in xmgrace:
Homology Modelling:
Run blast
blastp -query query.fasta -db pdbaa \
  -out blastp_results.tsv \
  -evalue 1e-5 \
  -outfmt "6 qcovs pident evalue sacc ssciname stitle" \
  -max_target_seqs 10 \
  -num_threads 4
output saved in blastp_results.tsv
chose  because it had highest %identity*querycoverage value.

Go to rcsb.pdb → download in legacy pdb format
kept only chain A
pdb_selchain -A 7VKQ.pdb > 7VKQ_A.pdb
opened in chimera to look for non standard residue
chimera → select → residue
mutated SAH back to Cystiene using pymol
wizard → muategenesis → protien → select residue → mutate to → apply → save pdb
deleted hetero atoms
pdbfixer input.pdb --add-atoms=all --add-residues=all --replace-nonstandard --output=fixed.pdb
pdb_delhetatm 7VKQ_mutated.pdb > 7VKQ_nohet.pdb
delete anisou record
grep -v '^ANISOU' 7VKQ_nohet.pdb > 7VKQ_noanisou.pdb
delete seqres just by backspace in text editor
delete alternate conformers
pdb_selaltloc -A 7VKQ_noanisou.pdb > 7VKQ_clean.pdb
Remove tail using chimera:
hover over start and end of the tail remove the tail residue from pdb file using backspace in file editor
save the sequence using chimera
favourite → sequence → file → save as sto
from the stockholm format take copy sequence in a new fasta file
perform pairwise sequence alignment
needle -asequence query.fasta -bsequence template.fasta   -gapopen 10 -gapextend 0.5   -outfile alignment.fasta   -aformat fasta -auto
create align.pir file
use python script got it from internet
make necessary changes
run modeller
mod10.7 script.py
after the run obtained 3 models
to find best model
chimera model1 template.pdb
tools → structure comaprison → matchmaker
repeat for all the 3 models

Generate ramachandran plot
launch ccp4
ccp4i → procheck → select model → uncheck sfcheck → run → run now
(/home/ibab/sem3/ST/pdb_manipulation)
based on RMSD values and Ramachandran plot the best predicted model is


Docking:

to find the residues we can use online servers like prankweb however there is some issue with my system because of which I am not able to upload files on web servers my file manager keeps on crashing so I tried and found another alernative fpocket web but here as well drag and drop feature did not work
installed fpocket and the command:
fpocket -f protein.pdb
 
chimera protien.pdb ligand.pdb
select → chain A → protein.pdb → action → surface → show
adt
file → read molecule → protein.pdb
edit → hydrogen → add → all hydrogens
grid → macrmolecule → choose → save as pdbqt
ligand → input → add ligand.pdb
added
ligand → ouput → save as → pdbqt
grid → set map type → choose ligand
grid → box →  cover entire protein
file → close → saving current (in grid options dialogue box only)
grid → output → save gpf
give file name grid.gpf
run → run autogrid →
write log file name
in the command line write the log file name before &
docking → macromolecule → set rigid file name → choose protein.pdbqt
docking → ligand →open → ligand.pdbqt → accept
docking → search parameters → genetic algorithm → accept
docking → docking parameters → accept
docking → output → lamarckian → lig1.dpf
run →run autodock →  dlg file will get selected → lig1.dlg should be log file name and give the log file name before the & as well → launch
right click on old ligand → hide
show protein as surface
analyze → docking → select log file
analyze → conformations → there are 10 different found
select the best
click on & → write complex → save
click on & → write current → save (this saves ligand)
obabel -ipdbqt complex.pdbqt -opdb -O pose1.pdb
make a config file for vina he did a google search
here this config.txt shall have parameters that are of the best docking site (after blind docking you found the best docking site after that make another grid around it write down the values AND THEN DO VINA WITH THOSE SO NOW YOU ARE DOING DOCKING AS PER THE SPECIFIC CORDINATES)
vina --ligand ligand.pdbqt --config config.txt --out run1.pdbqt
or run vina script
to visualize autodock result:
load complex.pdbqt → select ligand groups → define ligand ligand interactions → show 2D diagram
to visualize vina result:
adt → load protein.pdbqt →analyze → vina result → select the ligand.pdbqt → ok → using arrows you can browse different bindings

convert the run1.pdbqt to run1.pdb
obabel -ipdbqt run1.pdbqt -opdbqt -O run1.pdb
vina result have so many protien to take just one copy paste one model till endmdl in another pdbfile ligvina.pdb
launch discoverystudio
open receptor.pdb
file → insert → ligvina.pdb
receptorligand interactions → show interactions → click on check box for ligand → show 2D diagram






CCP4:
ccp4i
Click on popup after Refinement -> Select Coordinate Utilities
Edit PDB file
Job title : 2kho_clean_rm_alt
Use pdbcur to remove alternative conformations
PDB in : 2kho.pdb
PDB out : 2kho_pdbset1.pdb
Run Now
Edit PDB file
Job title : 2kho_clean_rm_aniso
Use pdbcur to remove anisotropic U's
PDB in : 2kho_pdbset1.pdb
PDB out : 2kho_pdbset2.pdb
Run Now
Ramachandran Plot:
Ramchandran plot analysis
ccp4i
Program list -> Procheck
Job title: cg rm-plot
Disable Run Sfcheck
Coords in: Full path : /home/ibab/NIYATI/SEM3/ST/15_9_25_energy_min/em-cg-mol.pdb
Run
View files from job -> cg rm-plot

Chimera:
select bond → tools → structure editing → build structure → adjust bonds → Delete
select atoms amongst which you want to make bond → tools → structure editing → build structure → adjust bonds → add
select atom → tools → structure editing → nuild structure → modify structure → atom name bond geometry
opened ligand2_acetyl.pdb → selected bonds about which rotations can be done(not part of ring) → tools → structure editing → adjust torsions → rotated torsion angles of all bonds to generate compact structure
file → select -> select all
tools -> structure analysis -> find clashes and contacts
yellow lines less than 2 A short contacts 
Favourites → CommandLog
Pymol:
sele chA, (chain A)
sele nag, (resn NAG)
sele atp_A, (chain A & resn ATP)

select molecule → surface → file → save → surface → file → open → A → generate → vacuum electrostatics → protien contact potential (local)

ball and sticks: s→ sphere → settings → edit all → sphere_scale 0.25

adjusting helix appearance:
settings → editall
cartoomn_sampling 10
cartoon_oval_width 0.15
cartoon_oval_length 0.8
cartoon_transparency 0.2
DynDom:
Extract chain A (pdb-tools)
pdb_selchain -A 4b9q.pdb > 4b9q_chainA.pdb
pdb_selchain -A 2kho.pdb > 2kho_chainA.pdb

Remove heteroatoms
pdb_delhetatm 4b9q_chainA.pdb > 4b9q_chainA_clean.pdb
pdb_delhetatm 2kho_chainA.pdb > 2kho_chainA_clean.pdb

Renumber residues from 1
pdb_reres -1 4b9q_chainA_clean.pdb > 4b9q_chainA_renum.pdb
pdb_reres -1 2kho_chainA_clean.pdb > 2kho_chainA_renum.pdb

Confirm max residue numbers (robust)
awk '/^ATOM/ {print substr($0,23,4)+0}' 4b9q_chainA_renum.pdb | sort -n | uniq | tail -n1
awk '/^ATOM/ {print substr($0,23,4)+0}' 2kho_chainA_renum.pdb | sort -n | uniq | tail -n1

Trim both to common region (example 1:585)
pdb_selres 1:585 4b9q_chainA_renum.pdb > 4b9q_final.pdb
pdb_selres 1:585 2kho_chainA_renum.pdb > 2kho_final.pdb

fix chain identifiers:
pdb_chain -A input.pdb > output.pdb
SSP:
Chou-Fasman, GOR and Neural network
Paste the query sequence
Tick all the check boxes
Click on start prediction

Discovery studio
File -> New -> Protein sequence window
Paste the query seq
Sequence -> Secondary structure -> Predict
If you have a pdb file then you can open that
Right click on the screen -> Sequence

JPred
Copy query sequence
Paste it into input sequence block
Give a fasta tag while running this
Remove dashes or extra space or new line character if any present
Click on Make Prediction
Will take onger duration so one should gve their email address

MD simulations:

VMD -> Visual Molecular Dynamics
- used to view trajectory file
1. Load the md.gro file
2. Load the md.xtc file (xtc is compressed and trr is uncompressed)
3. Selection: all -> Style: NewCartoon
4. Selection: protein -> Click on + -> Style: Lines

We need to get an idea about which two residues are interacting

Can change that to orthographic view (Display -> Orthographic)
Mouse -> Label -> Atoms

Analysis of MD Simulation results
1. Calculate RMSD
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -what -fit
On running the above command we got an Error

--> gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -what rmsd
2
2
We can compare entire protein / Protein without H / C-alpha / backbone
--> xmgrace rmsd.xvg
	In xmgrace
	Data -> Transformations -> Running average -> length of average: 10

2. Calculate RMSF (root mean square fluctuation)
To know whether there is binding site or not and which are the interacting residues

--> gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -oq bfac.pdb -res
	-oq -> rmsf values are converted to b factor  
	Select a group: 1
--> xmgrace rmsf.xvg
	In xmgrace: Plot -> Axis properties -> Start from 1

--> Open bfac.pdb file in pymol
	color -> spectrum -> b-factors

3. Calculate the radius of gyration
Used to quantify the compactness of protein
Measures how far the atoms of a molecule are from its center of mass

--> gmx gyrate  -f md.xtc -s md.tpr -o gyrate.xvg
	Select a group: 1
--> xmgrace gyrate.xvg

4. H-bonding
Choose 2 groups / residues and check whether H-bonding is happening or not
Need to create index file here as we need to specify two groups for analysis
Earlier we saw group list -> it is standard
Now we want to know whether there is H-bond between residue no 6 and 7, so we need to create new index

--> gmx make_ndx -f md.gro -o r6-r7.ndx
	r 6
	r 7
	q -> quit
--> gmx make_ndx -f md.gro -o r6-r7.ndx -n r6-r7.ndx
	17 r_6                 :    15 atoms
	18 r_7                 :    24 atoms
--> gmx hbond -f md.xtc -s md.tpr -num hbond.xvg -n r6-r7.ndx
	Select a group: 17
	Select a group: 18

5. Find the distance between two atoms
Calculates distance between pairs of selection
For eg. Find distance between C-alpha of 1st and 8th residue
For this we again need to make index -> A single group with two atoms

--> gmx make_ndx -f md.gro -o dist.ndx
	r 1 & a CA
	r 8 & a CA
	19 | 20 -> to merge two selections / groups into one group
	q
--> gmx distance -f md.xtc -s md.tpr -n dist.ndx -oall r1-r8-CA-dist.xvg
	Select a group: 21
	Ctrl D
xmgrace r1-r8-CA-dist.xvg

Ideally we should observe a stable plateau after equilibration
