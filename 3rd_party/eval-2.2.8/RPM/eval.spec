Summary: A System for Comparing and Analyzing Gene Sets
Name: eval
Version: 2.2.8
Release: 1
Group: Applications/System
Copyright: GPL
URL: http://genes.cse.wustl.edu/
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-root
Requires: perl
#BuildRequires: perl-libxml-enno
AutoReq: 0
BuildArch: noarch

%description
Eval is a Perl-based system for comparing sets of genes from any source,
such as pure computational predictions and manual annotations. Many types
of analyses can be run allowing for detailed, comprehensive comparisons. 
Results are displayed in a convenient, compact format.  

%package    common
Obsoletes: eval
Summary: A System for Comparing and Analyzing Gene Sets
Group: Applications/System

%package    X11
Summary: An X11 frontend to the Eval scripts
Group: Applications/System
Requires: eval-common perl-Tk

%description common
Eval is a Perl-based system for comparing sets of genes from any source,
such as pure computational predictions and manual annotations. Many types
of analyses can be run allowing for detailed, comprehensive comparisons. 
Results are displayed in a convenient, compact format.  

%description X11
An X11 frontend to the Eval scripts.

%prep
%setup -q 

%build

%install 
rm -rf $RPM_BUILD_ROOT
eval `perl '-V:installvendorlib'`
mkdir -p $RPM_BUILD_ROOT/$installvendorlib
mkdir -p $RPM_BUILD_ROOT/usr/bin

install -m0644 Eval.pm $RPM_BUILD_ROOT/$installvendorlib
install -m0644 GTF.pm $RPM_BUILD_ROOT/$installvendorlib
cp *.{pl,py} $RPM_BUILD_ROOT/usr/bin
find $RPM_BUILD_ROOT/usr -type f -not -name "eval.pl" -print | \
        sed "s@^$RPM_BUILD_ROOT@@g" > eval-common-%{version}-filelist 

mkdir -p $RPM_BUILD_ROOT/usr/lib/python2.3/site-packages
install -m0644 RangeFinder.py $RPM_BUILD_ROOT/usr/lib/python2.3/site-packages/RangeFinder.py

%clean
rm -rf $RPM_BUILD_ROOT

%post


%preun

%files common -f eval-common-%{version}-filelist
%defattr(-, root, root)
%doc help/*

%files X11 
%defattr(-, root, root)
/usr/bin/eval.pl
/usr/lib/python2.3/site-packages/RangeFinder.py

%changelog
* Tue Mar 29 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.8-1
- Adding chr22 example files

* Tue Mar 20 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.7-1
- Released version of eval 2.2.6 has a dependency on bioperl, bumping the 
    version to eliminate this dependency

* Tue Mar 20 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.6-1
- Fixed a bug where crash happens when not enough CDS is present for a full
    start codon.

* Fri Feb 23 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.5-1
- Added some support for GFF2 in gff3_to_gtf.pl, however
    there isn't much to recommend it.

* Tue Feb 20 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.4-1
- Corrected some bugs in the GTF.pm error correcting.

* Tue Feb 20 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.3-1
- Added methods to fix some gtf errors to GTF.pm in Transcript
- Added gff3_to_gtf.pl, which currently works for GMAP GFF3 (-f 2) output,
    provided that only one path per alignment is output. (Note this
    requires bioperl, but I don't want to make this a requirement for
    eval yet.)

* Mon Oct 25 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.2-1
- Bug fix to the create_utr_objects_from_exons subroutine

* Mon Jul 17 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2.1-1
- Made a fix in merge_gtf_transcripts.py to allow comments

* Mon Jul 17 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-10
- Adding the RangeFinder.py module

* Mon Jul 17 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-9
- Repackaging with RPM included.

* Mon Jul 17 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-8
- Moved obsoletes to eval-common only.

* Mon Jul 17 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-7
- Added the Obsoletes: bit.

* Sun Jul 16 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-6
- ?

* Fri Jun 02 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-5
- Split the package into two packages, X11 with eval.pl and common without.

* Fri Jun 02 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-5
- Made some slight changes to evaluate_gtf, eschewing some output
- Added filter_badlist.pl, since it's part of the annotation pipeline

* Fri Jun 02 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-4
- Fixed a bug in get_general_stats.pl where the ASEs were not marked if
    the input file is a list and not a gtf file.

* Wed May 10 2006 Bob Zimmermann <rpz@cse.wustl.edu> - 2.2-3
- Made some changes to the altsplicing detection to make it more efficient
- Added installation of the merge_gtf_transcripts.py script

* Wed May 10 2006 Brian Koebbe <koebbe@wustl.edu> - 2.2-2
- make it a noarch

* Tue May 9 2006 Brian Koebbe <koebbe@wustl.edu> - 2.2-1
- Initial build
