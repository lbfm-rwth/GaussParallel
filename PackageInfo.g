#
# GaussParallel: Parallel gaussian algorithm for finite fields
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "GaussParallel",
Subtitle := "Parallel gaussian algorithm for finite fields",
Version := "1.0.0dev",
Date := "05/05/2021", # dd/mm/yyyy format

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Jendrik",
    LastName := "Brachter",
    #WWWHome := "mailto:brachter@cs.uni-kl.de",
    Email := "brachter@cs.uni-kl.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Sergio",
    LastName := "Siccha",
    #WWWHome := "mailto:siccha@mathematik.uni-kl.de",
    Email := "siccha@mathematik.uni-kl.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := false,
    FirstNames := "Emma",
    LastName := "Ahrens",
    #WWWHome := "mailto:emma.ahrens@rwth-aachen.de",
    Email := "emma.ahrens@rwth-aachen.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
],

#SourceRepository := rec( Type := "TODO", URL := "URL" ),
#IssueTrackerURL := "TODO",
PackageWWWHome := "https://github.com/lbfm-rwth/GaussParallel/",
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL     := Concatenation( ~.PackageWWWHome,
                                 "archive/v", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  
  "The <span class=\"pkgname\">GaussParallel</span> package provides an \
  implementation of a parallel Gaussian elimination algorithm as described in \
  <a href=\"https://arxiv.org/abs/1806.04211\">\
  A parallel algorithm for Gaussian elemination over finite fields</a>.",

PackageDoc := rec(
  BookName  := ~.PackageName,
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := ~.Subtitle,
),

Dependencies := rec(
  GAP := ">= 4.10",
  # list of pairs [package name, version], package name is case
  # insensitive, exact version denoted with '=' prepended to version string.
  # without these, the package will not load
  NeededOtherPackages := [],
  SuggestedOtherPackages := [["io", "4.5.4"],
                             ["GAPDoc", "1.6.2"],
                             ["Gauss", "2018.09.08"],
                             ["AutoDoc", "2019.05.20"]],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


