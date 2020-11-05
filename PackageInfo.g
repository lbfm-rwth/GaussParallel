#
# GaussPar: Parallel gaussian algorithm for finite fields
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "GaussPar",
Subtitle := "Parallel gaussian algorithm for finite fields",
Version := "0.1",
Date := "24/10/2018", # dd/mm/yyyy format

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Jendrik",
    LastName := "Brachter",
    #WWWHome := TODO,
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
    #WWWHome := TODO,
    Email := "siccha@mathematik.uni-kl.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Emma",
    LastName := "Ahrens",
    #WWWHome := TODO,
    Email := "emma.ahrens@rwth-aachen.de",
    #PostalAddress := TODO,
    #Place := TODO,
    #Institution := TODO,
  ),
],

#SourceRepository := rec( Type := "TODO", URL := "URL" ),
#IssueTrackerURL := "TODO",
PackageWWWHome := "https://github.com/lbfm-rwth/GaussPar/",
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL     := Concatenation( ~.PackageWWWHome,
                                 "/", ~.PackageName, "-", ~.Version ),

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

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "GaussPar",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Parallel gaussian algorithm for finite fields",
),

Dependencies := rec(
  GAP := ">= 4.9",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


