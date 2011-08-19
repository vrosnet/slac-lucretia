/* file containing various stuff useful in the generation of
   messages during Lucretia lattice verification */

	char* zwf = "WF.ZSR" ;
	char* twf = "WF.TSR" ;

/* standard words commonly used: */

	char* KLYSStr = "KLYSTRON" ;
	char* KlysStr = "Klystron" ;
	char* PSStr = "PS" ;
	char* GIRDERStr = "GIRDER" ;
	char* GirderStr = "Girder" ;
	char* BEAMLINEStr = "BEAMLINE" ;
	char* ElemStr = "Element" ;
	char* ErrStr = "Error" ;
	char* WarnStr = "Warning" ;
	char* InfoStr = "Info" ;
	char* reqdStr = "required" ;
	char* optStr = "optional" ;
	char* TLRStr = "WF.TLR" ;
	char* TLRErrStr = "WF.TLRErr" ;

	char* Keywd ;
	char* KEYWD ;
	char* parname ;

/* standard messages */

/* missing parameter */

	char* MissPar = "%s: %s %d, %s parameter %s missing" ;

/* invalid length */

	char* BadLen = "%s: %s %d, parameter %s length incorrect" ;

/* points to non-existent table entry */

	char* NoSuch = "%s: %s %d points at non-existent %s %d" ;

/* inconsistent pointers */

	char* Inconsis = "%s: inconsistency between %s %d and %s %d" ;

/* zero elements */

	char* ZeroElt = "%s: %s %d points at zero elements" ;