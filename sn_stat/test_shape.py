import sn_stat as sn

def test_shapeana():
    B = sn.rate(1)
    S = sn.rate(2, range=[-1,1])
    det = sn.DetConfig(S=S,B=B)
    ana = sn.ShapeAnalysis([det])
    assert ana is not None
