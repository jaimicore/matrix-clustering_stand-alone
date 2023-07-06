#!/usr/bin/env python

'''NAME
        annotate-html-radialtree.py

VERSION
        %(version)s

AUTHOR
        Walter Santana-Garcia <walter.santana.garcia@ens.fr>

DESCRIPTION
        Creates a circos plot with a table annotation previously processed
        by annotate_matrix-clustering.R and a radial tree from RSAT::matrix-clustering.

CATEGORY
        matrix-tools

ARGUMENTS
    GENERAL OPTIONS
        --version             show program's version number and exit
        -h, --help            show this help message and exit
        -v #, --verbosity=#   set verbosity to level #


        -i #, --input=#       Mandatory option. The directory path result from
                              RSAT::matrix-clustering containing a radial tree
                              and an annotation processed previously.
'''

VERSION = '20190926'
USAGE = '''annotate-html-radialtree.py -i directory
                [-h | --help]
'''



########################################
#                                      #
# IMPORTS
#                                      #
########################################
import os
import re
import sys
import argparse
import time

########################################
#                                      #
# COMMAND LINE OPTIONS
#                                      #
########################################
parse = argparse.ArgumentParser()
parse.add_argument("-i", required=True, action="store",      dest="input")

args = parse.parse_args()


########################################
#                                      #
#   MAIN
#                                      #
########################################

# Rename variable
mclust_dir    = args.input
# Create prefix
mclust_prefix   = mclust_dir.strip("/").split("/")[-1]
mclust_dir      = os.path.dirname(mclust_dir)


# Strip absolute path
mclust_dir_tree = mclust_prefix + "_trees"

# Create file name for input hmtl
mclust_html_file = mclust_dir + "/" + mclust_prefix   + "_D3_radial_tree.html"
# Create file name for input annotation JSON
mclust_json_file = mclust_dir + "/" + mclust_dir_tree + "/annotation_matrix.json"
# Create file name for output html radial tree
mclust_out_file  = mclust_dir + "/" +  mclust_prefix  + "_tree_annotated.html"


# Test if matrix-clustering directory exists
if not os.path.exists(mclust_dir):
    print ("Unable to find matrix-clustring output directory {0}".format(mclust_dir))
# Test if matrix-clustering radial html tree exists
if not os.path.exists(mclust_html_file):
    print ("Unable to find file {0}".format(mclust_html_file))
# Test if matrix-clustering json annotation exists
if not os.path.exists(mclust_json_file):
    print ("Unable to find file {0}".format(mclust_json_file))


# Create regex expressions
regex_insert_variable   = re.compile(r"<script type=\"text/javascript\">\n$")
regex_save_textlabels   = re.compile(r"node\.append\(\"svg:a\"\)\n$")
regex_end_textlabels    = re.compile(r"^\s+;\n$")
regex_insert_annotation = re.compile(r"^\s+\.attr\(\"preserveAspectRatio\".+\n$")
#regex_insert_annotation = re.compile(r"^\}\)\;\n")

# Open filehandlers
#fh_annotation   = open(mclust_json_file, 'r')
fh_html         = open(mclust_html_file, 'r')
fh_output       = open(mclust_out_file,  'w')

# Load all the annotation file content in a variable
with open(mclust_json_file, 'r') as fh_annotation:
    annotation = fh_annotation.read().replace("\\","").replace("},","},\n")


# Strip 2 leading characters from annotation
annotation = annotation[2:]
# Strip 2 trailing characters from annotation
annotation = annotation[:-3]
# Add leading character  to annotation
annotation = "var data_sample = " + annotation
# Add trailing character  to annotation
annotation = annotation + ";\n"

# Create d3.js annotation functions
d3_code = ("\n\n// Add background elements for selection\n"
"     vis.selectAll('rect_selection')\n"
"        .data(data_sample)\n"
"        .enter()\n"
"        .append('path')\n"
"        .attr('class', 'rect1')\n"
"        .attr('id', function(d){ return( 'rect_sel' + d.id_motif) })\n"
"        .attr('d', d3.svg.arc()\n"
"          .startAngle(function(d){ return( (d.start * Math.PI)/180)  }  )  //converting from degs to radians\n"
"          .endAngle(function(d){   return( (d.end   * Math.PI)/180)  }  ) //just radians\n"
"          .innerRadius(innerRad_start)         // This is the size of the donut hole\n"
"          .outerRadius(innerRad_end)\n"
"        )\n"
"        .style('fill',function(d){return(d.colour)})\n"
"        .attr('transform', 'translate(0,0)')\n"
"        .style('stroke-width', '1px')\n"
"        .style('opacity', 0.15);\n"
"\n"
"            \n"
"        // Color algorithm Non-Validated text layer\n"
"            vis.selectAll('annotations')\n"
"                .data(data_sample)\n"
"                .enter()\n"
"                .append('path')\n"
"                .attr('class','annotation1')\n"
"		  .attr('id', function(d,i){return('path' + i) })\n"
"                .attr('d', d3.svg.arc()\n"
"                  .startAngle(function(d){ return( d.start * (Math.PI/180) ) }  )  //converting from degs to radians\n"
"                  .endAngle(function(d){ return( d.end * (Math.PI/180) ) }  ) //just radians\n"
"                  .innerRadius(innerRad_end)         // This is the size of the donut hole\n"
"                  .outerRadius(innerRad_end + (30*1))\n"
"                )\n"
"                .attr('fill', function(d){return(d.colour)})\n"
"                .attr('stroke', 'white')\n"
"                .attr('transform', 'translate(0,0)')\n"
"                .style('stroke-width', '2px')\n"
"                .style('opacity', 1);\n"
"\n"
"\n"
"        // Color algorithm Non-Validated text layer\n"
"        vis.selectAll('annotation1')\n"
"            .data(data_sample)\n"
"            .enter()\n"
"            .append('text')\n"
"            .attr('dy', 20)\n"
"            .attr('x',33)\n"
"            .style('font-size', '15px')\n"
"            .append('textPath')\n"
"            .attr('startOffset','50%')\n"
"            .style('text-anchor','middle')\n"
"            //.attr('stroke','black')\n"
"            //.attr('fill','black')\n"
"            .attr('xlink:href', function(d,i){return('#path' + i) })\n"
"            .text(function(d){\n"
"              var token = d.motif_id.split('_').slice(-2).slice(0);\n"
"              token = token[0];\n"
"              if(/^UN/.test(token)){\n"
"                var text_content = d.class_nb + '*';\n"
"              } else {\n"
"                var text_content = d.class_nb;\n"
"              }\n"
"              return(text_content);\n"
"            })\n"
"            .attr('font-family', 'Ubuntu Mono')\n"
"            ;\n"
"\n"
"\n"
"            vis.selectAll('g.text')\n"
"                .data(nodes)\n"
"                .enter()\n"
"                .append('g')\n"
"                .attr('class', 'node_text')\n"
"                .attr('transform', function(d) { return 'rotate(' + (d.x - 90) + ')translate(' + d.y + ')'; })\n"
"                .append('svg:a')\n"
"//                .attr('xlink:href', function(d) { return d.link_ext; })\n"
"                .attr('xlink:href', 'https://jaspar.uio.no/matrix/MA0261.1/')\n"
"                .attr('target', '_blank')\n"
"                .append('text')\n"
"                .text(function(d) { return d.children ? '' :  d.name; })\n"
"                .attr('fill', function(d) { return d.color; })\n"
"                .style('font-size', '12px')\n"
"                .attr('font-family', 'Ubuntu Mono')\n"
"                .attr('dx', function(d) { return d.x < 180 ? 10 : -10; })\n"
"                .attr('dy', '.31em')\n"
"                .attr('text-anchor', function(d) { return d.x < 180 ? 'start' : 'end'; })\n"
"                .attr('transform', function(d) { return d.x < 180 ? null : 'rotate(180)'; });\n")


#  return d.color;  This color must be the class color

# Create flag
save = 0
# Create empty string
chunksave = ""

# Skip comment lines
for line in fh_html:
    # Write json as variable to new file
    if regex_insert_variable.search(line):
        fh_output.write(line)
        # Insert variable to file
        fh_output.write(annotation)
    # Insert code
    elif regex_insert_annotation.search(line):
        fh_output.write(line)
        # Insert d3.js code functions
        fh_output.write(d3_code)
    # Start saving chunk
    elif regex_save_textlabels.search(line):
        # Create flag for chunk saving
        save = 1
        # Start saving chunk
        chunksave =  line
    # End saving chunk
    elif regex_end_textlabels.search(line):
        # End chunk saving
        save = 0
        # Add last line to chunk
        chunksave = chunksave + line
    # Save chunk line
    elif save:
        # Keep saving chunk
        chunksave = chunksave + line
    else:
        fh_output.write(line)

# Close filehandlers
fh_html.close()
fh_output.close()
# Verbose end
print("; Radial tree html annotated succesfully")
