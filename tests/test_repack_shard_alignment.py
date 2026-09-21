"""Regression: edge dimensions must still produce valid whole-chunk shards."""
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import zarr


@pytest.mark.parametrize('shape,chunks',[( (3,2),(64,1000)), ((193,5),(64,2))])
def test_repack_preserves_edge_arrays(tmp_path,shape,chunks):
    src=tmp_path/'input.zarr';dst=tmp_path/'output.zarr'
    group=zarr.open_group(str(src),mode='w',zarr_format=3)
    values=np.arange(np.prod(shape),dtype=np.int16).reshape(shape)
    group.create_array('call_DP',data=values,chunks=chunks)
    group.create_array('variant_position',data=np.arange(shape[0],dtype=np.int32),chunks=(chunks[0],))
    script=Path(__file__).parents[1]/'scripts/repack_vcz_zarr_v3.py'
    result=subprocess.run([sys.executable,str(script),str(src),str(dst),'--workers','1',
                    '--variant-shard-chunks','25','--compressed-passthrough'],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    output=zarr.open_group(str(dst),mode='r')['call_DP']
    np.testing.assert_array_equal(output[:],values)
    assert output.chunks==chunks
    assert all(s%c==0 for s,c in zip(output.shards,output.chunks))
