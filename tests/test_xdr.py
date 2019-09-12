#!/usr/bin/env python
import pytest
from pmx.model import Model
from pmx.xtc import Trajectory
from numpy.testing import assert_almost_equal

@pytest.mark.parametrize("struct, traj", [
    ('peptide.pdb', 'peptide.trr'),
    ('peptide.pdb', 'peptide.xtc')
])
def test_trajectory_xdr(gf, struct, traj):
    m = Model(gf(struct))
    t = Trajectory(gf(traj))
    for frame in t:
        frame.update( m )
    assert_almost_equal(m.atoms[2].x[1], 34.70, decimal=2)
    t.close()
    
