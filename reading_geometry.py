def get_geometry(domain):
    """ Reads a geometry from a file.
    Args:
      domain(str) : type of domain to open."""
    if (domain == 0) :
        ls_domain = "../domains/2squares_inv.xml"#n64x64_p3x3.xml"
        import caid.cad_geometry as cg
        geo = cg.cad_geometry(ls_domain)
    elif (domain == 1) :
        ls_domain = "../domains/disk_100100.xml"
        import caid.cad_geometry as cg
        geo = cg.cad_geometry(ls_domain)
        # from caid.cad_geometry import circle
        # geo = circle()
    elif (domain == 2) :
        ls_domain = "../domains/disk_5mp_n64x64_p3x3.xml"
        import caid.cad_geometry as cg
        geo = cg.cad_geometry(ls_domain)
    elif (domain == 3) :
        ls_domain = "../domains/pinched_Disk-5mp_n32x32_p3x3.xml"
        import caid.cad_geometry as cg
        geo = cg.cad_geometry(ls_domain)
    elif (domain == 4) :
        ls_domain = "../domains/4squares.xml"
        import caid.cad_geometry as cg
        geo = cg.cad_geometry(ls_domain)
    elif (domain == 5) :
        ls_domain = "../domains/annulus.xml"
        import caid.cad_geometry as cg
        geo = cg.cad_geometry(ls_domain)
    else:
        raise SystemExit("ERROR in globals_variables:" \
                         + " no domain associated to that number")
    return geo

def get_patches(geo):
    list_p = []
    list_n = []
    for nrb in geo :
        list_p.append(nrb.degree)
        list_n.append(nrb.shape)

    print "list p:", list_p
    print "list n:", list_n
    p = list_p[0][0]
    n = list_n[0][0]
    print(" n = "+str(n)+" \n p = "+str(p))

    npatchs = geo.npatchs
    list_patchs = list(range(0, npatchs))
    print(" Num of patchs : " + str(npatchs))
    print("")

    lpi_ordregl = (p,p)
    return list_patchs
