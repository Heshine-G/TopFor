import subprocess
import os
from modules.remove import process_mol2_file

def run_antechamber_for_all(mol2_files):
    for input_mol2 in mol2_files:
        base_name = os.path.splitext(os.path.basename(input_mol2))[0]
        residue_name = base_name.upper()
        output_dir = base_name
        os.makedirs(output_dir, exist_ok=True)

        ac_output = os.path.join(output_dir, f"{base_name}.ac")
        mol2_output = os.path.join(output_dir, f"{base_name}.mol2")
        lib_output = os.path.join(output_dir, f"{base_name}.lib")
        prepin_output = os.path.join(output_dir, f"{base_name}.prepin")
        mc_output = os.path.join(output_dir, f"{base_name}.mc")
        frcmod_output = os.path.join(output_dir, f"{base_name}.frcmod")
        gaff_frcmod_output = os.path.join(output_dir, f"{residue_name}_gaff.frcmod")
        ff14SB_frcmod_output = os.path.join(output_dir, f"{residue_name}_ff14SB.frcmod")

        try:
            subprocess.run(f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo ac -o {ac_output} -c gas -at amber", shell=True, check=True)
            subprocess.run(f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo mol2 -o {mol2_output} -c gas -at amber", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Antechamber failed: {e}")
            continue

        leap_script = f"""
source leaprc.gaff
{residue_name} = loadmol2 {mol2_output}
edit {residue_name}
desc {residue_name}
set {residue_name} head {residue_name}.1.N
set {residue_name} tail {residue_name}.1.C
saveoff {residue_name} {lib_output}
quit
"""
        leap_file = os.path.join(output_dir, "leap.in")
        with open(leap_file, 'w') as f:
            f.write(leap_script)
        subprocess.run(f"tleap -f {leap_file}", shell=True, check=True)

        try:
            process_mol2_file(mol2_output, mc_output)
        except Exception as e:
            print(f"MC generation failed: {e}")
            continue

        subprocess.run(f"prepgen -i {ac_output} -o {prepin_output} -m {mc_output} -rn {residue_name}", shell=True, check=True)
        subprocess.run(f"parmchk2 -i {prepin_output} -f prepi -o {frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat", shell=True, check=True)
        subprocess.run(f"parmchk2 -i {ac_output} -f ac -o {gaff_frcmod_output}", shell=True, check=True)
        subprocess.run(f"parmchk2 -i {ac_output} -f ac -o {ff14SB_frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat", shell=True, check=True)

