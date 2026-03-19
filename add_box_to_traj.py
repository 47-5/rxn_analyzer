from ase.io import read, write


def add_box_to_traj(traj_filepath, ref_cell=None, ref_cell_filepath=None, new_traj_filepath='new.extxyz'):

    trajs = read(traj_filepath, index=':')

    if ref_cell is None:
        ref_cell = read(ref_cell_filepath).cell

    for atoms in trajs:
        atoms.cell = ref_cell

    write(new_traj_filepath, trajs)
    return trajs


if __name__ == '__main__':

    add_box_to_traj(
        traj_filepath='GaOH2__1_DMF_1_C2H4_rank01-pos-1.xyz',
        ref_cell_filepath='Ga__1_DMF_1_C2H4_rank02.cif',
        new_traj_filepath='GaOH2__1_DMF_1_C2H4_rank01-pos-1.extxyz'
    )


