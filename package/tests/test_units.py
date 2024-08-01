import sails

def test_pair():
    p = sails.Pair(1, 2)
    assert p.value1 == 1
    assert p.value2 == 2
    p.value1 = 3
    assert p.value1 == 3

def test_hkl():
    r = sails.HKL(1, 2, 3)
    assert r.h == 1
    assert r.k == 2
    assert r.l == 3
    r.h = 10
    assert r.h == 10

def test_reflection():
    i = sails.HKL(1, 2, 3)
    p = sails.Pair(1e-3, 1e8)

    r = sails.Reflection(i, p)
    assert r.hkl.h == i.h
    assert r.hkl.k == i.k
    assert r.hkl.l == i.l

    assert r.f_sigf.value1 == p.value1
    assert r.f_sigf.value2 == p.value2

    fwt = sails.Pair(1, 4)
    delfwt = sails.Pair(2, 5)

    r = sails.Reflection(i, p, fwt, delfwt)
    assert r.hkl.h == i.h
    assert r.hkl.k == i.k
    assert r.hkl.l == i.l

    assert r.f_sigf.value1 == p.value1
    assert r.f_sigf.value2 == p.value2
    assert r.fwt_phwt.value1 == fwt.value1
    assert r.fwt_phwt.value2 == fwt.value2
    assert r.delfwt_phdelwt.value1 == delfwt.value1
    assert r.delfwt_phdelwt.value2 == delfwt.value2

