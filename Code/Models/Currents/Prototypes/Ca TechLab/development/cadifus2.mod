TITLE Calcium ion accumulation and diffusion
: The internal coordinate system is set up in PROCEDURE coord_cadifus()
: and must be executed before computing the concentrations.
: The scale factors set up in this procedure do not have to be recomputed
: when diam or DFree are changed.
: The amount of calcium in an annulus is ca[i]*diam^2*vol[i] with
: ca[0] being the second order correct concentration at the exact edge
: and ca[NANN-1] being the concentration at the exact center

NEURON {
	SUFFIX cadifus2
	USEION ca READ cao, cai, ica WRITE cai, ica
	GLOBAL vol, Buffer0
	RANGE ipump
}
DEFINE NANN  4

UNITS {
        (mol)   = (1)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	FARADAY = (faraday)	 (10000 coulomb)
	PI	= (pi) (1)
}

PARAMETER {
	DFree = .6	(um2/ms)
	diam		(um)
	cao		(mM)
	ica		(mA/cm2)
	k1buf = 50	(/mM-ms)
	k2buf = 5	(/ms)
        k1=1.e10            (um3/s)
        k2=50.e7            (/s)	: k1*50.e-3
        k3=1.e10            (/s)	: k1
        k4=5.e6	            (um3/s)	: k1*5.e-4
	area		(um2)
} 
CONSTANT { volo=1  (liter)}

ASSIGNED {
	cai		(mM)
	vol[NANN]	(1)	: gets extra cm2 when multiplied by diam^2
	ipump           (mA/cm2)
	last_ipump           (mA/cm2)

}

STATE {
	ca[NANN]	(mM) <1.e-5> : ca[0] is equivalent to cai
	CaBuffer[NANN]  (mM)
	Buffer[NANN]    (mM)
        pump            (mol/cm2) <1.e-15>
        pumpca          (mol/cm2) <1.e-15>

}

LOCAL totpump, kd,totbuf

INITIAL {
           totpump=0.2
           pump=totpump/(1+1.e-18*k4*cao/k3)
           pumpca=2.e-22
	   ipump=0

           totbuf=1.2
           kd=k2buf/k1buf
           FROM i=0 TO NANN-1 {
                ca[i] = cai
		CaBuffer[i] =(totbuf*ca[i])/(kd+ca[i])
		Buffer[i] = totbuf - CaBuffer[i]
                }

}

BREAKPOINT {
	SOLVE state METHOD sparse
	last_ipump=ipump
	ica = ipump
}

LOCAL coord_done

INITIAL {
	if (coord_done == 0) {
		coord_done = 1
		coord()
	}
	: note Buffer gets set to Buffer0 automatically
	: and CaBuffer gets set to 0 (Default value of CaBuffer0) as well
	FROM i=0 TO NANN-1 {
		ca[i] = cai
	}
}

LOCAL frat[NANN] 	: gets extra cm when multiplied by diam

PROCEDURE coord() {
	LOCAL r, dr2
	: cylindrical coordinate system  with constant annuli thickness to
	: center of cell. Note however that the first annulus is half thickness
	: so that the concentration is second order correct spatially at
	: the membrane or exact edge of the cell.
	: note ca[0] is at edge of cell
	:      ca[NANN-1] is at center of cell
	r = 1/2					:starts at edge (half diam)
	dr2 = r/(NANN-1)/2			:half thickness of annulus
	vol[0] = 0
	frat[0] = 2*r
	FROM i=0 TO NANN-2 {
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2	:interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)	:exterior edge of annulus
					: divided by distance between centers
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2	:outer half of annulus
	}
}

LOCAL dsq, dsqvol : can't define local variable in KINETIC block or use
		:  in COMPARTMENT
KINETIC state {
	COMPARTMENT i, diam*diam*vol[i]*1(um) {ca CaBuffer Buffer}
        COMPARTMENT (1.e10)*area {pump pumpca}
        COMPARTMENT (1.e15)*volo {cao}

	~ ca[0] << (-(ica-last_ipump)*PI*diam*frat[0]*1(um)/(2*FARADAY))
	FROM i=0 TO NANN-2 {
		~ ca[i] <-> ca[i+1] (DFree*frat[i+1]*1(um), DFree*frat[i+1]*1(um))
	}
	dsq = diam*diam*1(um)
	FROM i=0 TO NANN-1 {
		dsqvol = dsq*vol[i]
		~ ca[i] + Buffer[i] <-> CaBuffer[i] (k1buf*dsqvol,k2buf*dsqvol)
	}
        ~ca[0] + pump <-> pumpca ((1.e-11)*k1*area, (1.e7)*k2*area)
        ~pumpca       <-> pump + cao ((1.e7)*k3*area, (1.e-11)*k4*area)

        ipump = 2*FARADAY*(f_flux-b_flux)/area

	cai = ca[0]
}
	
COMMENT
At this time, conductances (and channel states and currents are
calculated at the midpoint of a dt interval.  Membrane potential and
concentrations are calculated at the edges of a dt interval.  With
secondorder=2 everything turns out to be second order correct.
ENDCOMMENT

