#!/usr/bin/perl -w
#
# Generate an html contents file for the current version of fjcontrib
#
# To be run from the main directory of a checkout of fjcontrib.
# Best run from a directory that is a clean checkout of a tag.
#
# set to 1 to sort contribs alphabetically, 0 otherwise
$sort=1;
# set to 1 to include release date taken from svn tags, 0 otherwise.
$include_date=1;
# set to 1 to include dependencies
$include_deps=1;
# set to 1 to include minimal fastjet version
$include_minFJ=1;
# set to 1 to generate the page for LPTHE rather than for HEPFORGE
$lpthe=0;

$versions="contribs.svn";

# get svn info from one common location
if (! exists $ENV{'svn_read'}) {
  print STDERR '$svn_read environment variable must be set, but is not (source common.sh first?)',"\n";
  exit(-1)
} else {
  $svn_read=$ENV{'svn_read'};
  print STDERR '$svn_read is: ', $svn_read,"\n";
}

$svn=$svn_read."/contribs/";
# link to browse the svn
$svnBrowse='https://phab.hepforge.org/source/fastjetsvn/browse/contrib/contribs/';
# this ensures that one doesn't get the whole "blame" info, which is ugly
$svnPost='?as=source&blame=off';

$topversion=`head -1 VERSION`;
chomp $topversion;

# read in contribs.svn file, fill contribs hash
open (VERSIONS, "<$versions") || die "Could not open $versions";
%contribs_hash=();
$contribs_array=();
while ($line = <VERSIONS>) {
  if ($line =~ /^\s*([a-z][^\s]*)\s+([^\s]*)/i) {
    $contrib = $1;
    $version = $2;
    push @contribs_array, $contrib;
    $contribs_hash{$contrib} = $version;
  }
}

# sort contribs by alphabetical order
if ($sort) { @contribs_array = sort keys %contribs_hash; }

# write out html table
$list='';
foreach ( @contribs_array ) { 
    $contrib = $_;
    $version = $contribs_hash{$_};
    if ($version =~ /^[0-9]/) {$version = "tags/$version";}
    ($textversion = $version) =~ s/tags\///;
    if($include_date) {
      # extract date of last version tag from svn
      print STDERR "getting date for $svn$contrib/$version\n";
      $date = `svn info $svn$contrib/$version | grep "Last Changed Date" | awk '{print \$4}'`;
      # One could  also use the --xml option and parse appropriately the output:
      # $date = `svn --xml list $svn$contrib/$version`;
      # At this stage, XML parsing is not implemented though.
      #
      #print $contrib." ".$date."\n";
    }
    # extract dependencies and minimal fastjet version from config file
    $deps=""; $minFJ="";
    if ( -e "$contrib/FJCONTRIB.cfg" ) {
       $deps = `grep 'dependencies:' $contrib/FJCONTRIB.cfg`;
       $minFJ = `grep 'minimal_fastjet_version:' $contrib/FJCONTRIB.cfg`;
       # deps
       if ($deps) {
         chomp $deps;
         $deps =~ /^\s*dependencies:\s*(.*)/;
         $deps = $1;
	 if (! $deps ) { $deps = ""; } # to avoid a warning
       } 
       # minFJ
       if ($minFJ) {
         chomp $minFJ;
         $minFJ =~ /^\s*minimal_fastjet_version:\s*(.*)/;
         $minFJ = $1;
	 if (! $minFJ ) { $minFJ = ""; } # to avoid a warning
       }
    }   
    $list .= "<tr> <td class=\"contribname\"> 
                   <a href=\"$svnBrowse$contrib/$version/\">$contrib</a>
               </td> <td style=\"{text-align:center;}\"> $textversion </td>";
    if ($include_date) {$list .= "<td>$date</td>";}
    $list .=  "<td>";
    if (-e "$contrib/README") {
      $list .= '<a href="'.$svnBrowse.$contrib.'/'.$version.'/README'.$svnPost.'">README</a> ';
    }
    if (-e "$contrib/NEWS") {
      $list .= '<a href="'.$svnBrowse.$contrib.'/'.$version.'/NEWS'.$svnPost.'">NEWS</a> ';
    }
    if ($include_deps) {$list .= "<td>$deps</td>";}
    if ($include_minFJ) {$list .= "<td>$minFJ</td>";}
    $list .= "</td></tr>\n";
}

# create head
$head='
<html>
<head>
';

if ( $lpthe ) {
$head .= '<link rel="stylesheet" href="/style.css" type="text/css" media="screen"/>
<!--#include virtual="/commonheader.html"-->
'; }

$head .= '<style type="text/css">
.contriblist table,.contriblist tr,.contriblist td,.contriblist th {
  text-align:center;
  border:1px solid white;
  padding:4px;
  padding-right:6px;
}
td.contribname {
  background-color:#eeeeee;
  text-align:left;
}
th.contriblist {
  text-align:center;
  padding:6px;
  padding-left:10px;
  padding-right:10px;
  background-color:#dddddd;
}
td.spanned {
  text-align:center;
  background-color:#dddddd;
}
.triangles { 
  display: inline-block; 
  line-height: 1em; 
  text-align: center; 
  font-size: 0.5em; /* Adjust this to make the triangles smaller or larger */} 
.triangles span { display: block; } 
</style>
<script>
function sortTable(n) {
  var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
  table = document.getElementById("fjcontribtable");
  switching = true;
  //Set the sorting direction to ascending:
  dir = "asc"; 
  /*Make a loop that will continue until
  no switching has been done:*/
  while (switching) {
    //start by saying: no switching is done:
    switching = false;
    rows = table.rows;
    /*Loop through all table rows (except the
    first, which contains table headers):*/
    for (i = 1; i < (rows.length - 1); i++) {
      //start by saying there should be no switching:
      shouldSwitch = false;
      /*Get the two elements you want to compare,
      one from current row and one from the next:*/
      x = rows[i].getElementsByTagName("TD")[n];
      y = rows[i + 1].getElementsByTagName("TD")[n];
      /*check if the two rows should switch place,
      based on the direction, asc or desc:*/
      if (dir == "asc") {
         if (x.innerHTML > y.innerHTML) {
          //if so, mark as a switch and break the loop:
          shouldSwitch= true;
          break;
        } 
      } else if (dir == "desc") {
        if (x.innerHTML < y.innerHTML) {
          //if so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
     }
    }
    if (shouldSwitch) {
      /*If a switch has been marked, make the switch
      and mark that a switch has been done:*/
      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
      switching = true;
      //Each time a switch is done, increase this count by 1:
      switchcount ++;      
    } else {
      /*If no switching has been done AND the direction is "asc",
      set the direction to "desc" and run the while loop again.*/
      if (switchcount == 0 && dir == "asc") {
        dir = "desc";
        switching = true;
      }
    }
  }
}
</script>
</head>
<body>
';

if ( $lpthe ) {
$head .= '<!--#include virtual="/fjversion.html"-->
<div id="page">
    <div id="mast">
      <ul id="mainMenu">
        <li class=""><a href="/">Home</a></li>
	<li class=""><a href="/about.html">About</a></li>
 	<li class=""><a href="/all-releases.html">Releases</a></li>
        <li class=""><a href="/quickstart.html" >Quick start</a></li>
        <li class=""><a href="/repo/fastjet-doc-<!--#echo var="fjversion"-->.pdf" >Manual</a></li>
        <li class=""><a href="/repo/doxygen-<!--#echo var="fjversion"-->/" >Doxygen</a></li>
        <li class=""><a href="/tools.html" >Tools</a></li>
        <li class="selected"><a href="/contrib/" >Contrib</a></li>
 	<li class=""><a href="/FAQ.html" >FAQ</a></li>
      </ul>
    </div>
<div id="header">
</div>
<p><br>
'; }
 
$head .= 'Version '.$topversion.' of FastJet Contrib is distributed with the following packages<p>
';

if ($lpthe ) { 
$head .= '<p><br>
'; }

$head.='<table class="contriblist" id="fjcontribtable">
<tr><th class="contriblist" onclick="sortTable(0)">Package &nbsp; <div class="triangles"> <span>&#9650;</span> <span>&#9660;</span> </div></th> 
    <th class="contriblist">Version</th>';
if($include_date) {$head .= '
    <th class="contriblist" onclick="sortTable(2)">Release date &nbsp; <div class="triangles"> <span>&#9650;</span> <span>&#9660;</span> </div></th>';}
$head .= '
    <th class="contriblist">Information</th>';
if ($include_deps) {$head .= '
    <th class="contriblist">Dependencies</th>';}   
if ($include_minFJ) {$head .= '
    <th class="contriblist">Minimal FJ version</th>';}   
$head .= '</tr> 
';

# create tail
$tail='
</table>
';

if ( $lpthe) {
$tail .= '<!--#include virtual="/footer.html"-->
</div> <!-- page -->
'; }

$tail .= '</body>
</html>
';

# put everything together
print $head, $list, $tail;
