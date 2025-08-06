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

        leap_script = f"""
source leaprc.gaff
{residue_name} = loadmol2 {mol2_output}
set {residue_name} head {residue_name}.1.N
set {residue_name} tail {residue_name}.1.C
saveoff {residue_name} {lib_output}
quit
"""
        leap_file = os.path.join(output_dir, "leap.in")

        def run_cmd(cmd, error_msg):
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"{error_msg}: {result.stderr.strip()}")
                return False
            return True

        # 1) Generate AC
        ac_cmd = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo ac -o {ac_output} -c bcc -at amber"
        if not run_cmd(ac_cmd, "Antechamber AC failed"): continue
        print(f"Successfully generated AC file: {ac_output}")

        # 2) Generate MOL2
        mol2_cmd = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo mol2 -o {mol2_output} -c bcc -at amber"
        if not run_cmd(mol2_cmd, "Antechamber MOL2 failed"): continue
        print(f"Successfully generated MOL2 file: {mol2_output}")

        # 3) Run tleap
        with open(leap_file, 'w') as f:
            f.write(leap_script)
        if not run_cmd(f"tleap -f {leap_file}", "tleap failed"): continue
        print(f"Successfully ran tleap and created library: {lib_output}")

        # 4) MC generation 
        try:
            process_mol2_file(mol2_output, mc_output)
            print(f"Successfully generated MC file: {mc_output}")
        except Exception as e:
            print(f"MC generation failed: {e}")
            continue

        # 5) Prepgen
        prep_cmd = f"prepgen -i {ac_output} -o {prepin_output} -m {mc_output} -rn {residue_name}"
        if not run_cmd(prep_cmd, "prepgen failed"): continue
        print(f"Successfully generated prepi file: {prepin_output}")

        # 6) parmchk2 steps
        if run_cmd(f"parmchk2 -i {prepin_output} -f prepi -o {frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat", "parmchk2 prepi failed"):
            print(f"Successfully generated FRCMOD: {frcmod_output}")
        else:
            continue

        if run_cmd(f"parmchk2 -i {ac_output} -f ac -o {gaff_frcmod_output}", "parmchk2 gaff failed"):
            print(f"Successfully generated GAFF FRCMOD: {gaff_frcmod_output}")
        else:
            continue

        if run_cmd(f"parmchk2 -i {ac_output} -f ac -o {ff14SB_frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/parm10.dat", "parmchk2 ff14SB failed"):
            print(f"Successfully generated FF14SB FRCMOD: {ff14SB_frcmod_output}")
        else:
            continue

        
        print(f"\n\033[1mParameter Generation Successful for: {base_name}\033[0m")

        
