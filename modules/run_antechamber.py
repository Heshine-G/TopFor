import subprocess
import os
from modules.remove import process_mol2_file

def run_antechamber_for_all(mol2_files, backbone='ff19SB', sidechain='gaff2', charge='bcc'):
    backbone_parm_map = {
        'ff14SB': 'parm10.dat',
        'ff19SB': 'parm19.dat',
        'ff99SB': 'parm99.dat'
    }

    for input_mol2 in mol2_files:
        base_name = os.path.splitext(os.path.basename(input_mol2))[0]
        # Extract only the residue name (remove chain ID, position)
        residue_name = base_name.split('_')[0].upper()
        output_dir = residue_name
        os.makedirs(output_dir, exist_ok=True)

        ac_output = os.path.join(output_dir, f"{residue_name}.ac")
        mol2_output = os.path.join(output_dir, f"{residue_name}.mol2")
        lib_output = os.path.join(output_dir, f"{residue_name}.lib")
        prepin_output = os.path.join(output_dir, f"{residue_name}.prepin")
        mc_output = os.path.join(output_dir, f"{residue_name}.mc")
        frcmod_output = os.path.join(output_dir, f"{residue_name}.frcmod")
        gaff_frcmod_output = os.path.join(output_dir, f"{residue_name}_{sidechain}.frcmod")
        backbone_frcmod_output = os.path.join(output_dir, f"{residue_name}_{backbone}.frcmod")

        leap_script = f"""
source leaprc.{sidechain}
source leaprc.protein.{backbone}
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

        # Use dynamic charge model in both antechamber steps
        ac_cmd = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo ac -o {ac_output} -c {charge} -at amber"
        if not run_cmd(ac_cmd, "Antechamber AC failed"): continue

        mol2_cmd = f"antechamber -fi mol2 -i {input_mol2} -bk {residue_name} -fo mol2 -o {mol2_output} -c {charge} -at amber"
        if not run_cmd(mol2_cmd, "Antechamber MOL2 failed"): continue

        with open(leap_file, 'w') as f:
            f.write(leap_script)
        if not run_cmd(f"tleap -f {leap_file}", "tleap failed"): continue

        try:
            process_mol2_file(mol2_output, mc_output)
        except Exception as e:
            print(f"MC generation failed: {e}")
            continue

        prep_cmd = f"prepgen -i {ac_output} -o {prepin_output} -m {mc_output} -rn {residue_name}"
        if not run_cmd(prep_cmd, "prepgen failed"): continue

        parm_file = backbone_parm_map.get(backbone, 'parm10.dat')
        parmchk_backbone = f"parmchk2 -i {ac_output} -f ac -o {backbone_frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/{parm_file}"
        parmchk_sidechain = f"parmchk2 -i {ac_output} -f ac -o {gaff_frcmod_output} -a Y -p $AMBERHOME/dat/leap/parm/{sidechain}.dat"

        if run_cmd(parmchk_backbone, "parmchk2 backbone failed"):
            print(f"Successfully generated backbone FRCMOD: {backbone_frcmod_output}")
        else:
            continue

        if run_cmd(parmchk_sidechain, "parmchk2 sidechain failed"):
            print(f"Successfully generated sidechain FRCMOD: {gaff_frcmod_output}")
        else:
            continue

        print(f"\n\033[1mParameter Generation Successful for: {residue_name}\033[0m")
