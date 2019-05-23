/* $Id: filesuffix.h 213998 2018-03-03 04:59:34Z twu $ */
#ifndef FILESUFFIX_INCLUDED
#define FILESUFFIX_INCLUDED


#ifdef PMAP

#define FWD_FILESUFFIX "pf"
#define REV_FILESUFFIX "pr"

#else
#define IDX_FILESUFFIX "ref"
#endif

#define OFFSETS_FILESUFFIX "offsets"
#define POSITIONS_HIGH_FILESUFFIX "positionsh"
#define POSITIONS_FILESUFFIX "positions"

#define LOCAL_OFFSETS_FILESUFFIX "locoffsets"
/* #define LOCAL_POSITIONS_HIGH_FILESUFFIX "locpositionsh" */
#define LOCAL_POSITIONS_FILESUFFIX "locpositions"

#endif

