
********************************************************************************************************************************************************************
See COPYRIGHT.TXT for the copyright notice.

This folder contains an implementation of the algorithm to find a "maxmin-fair" decomposition for bipartite matching as described in: 
the paper "Fair-by-design matching" by David Garcia-Soriano and Francesco Bonchi.

Equivalently, it contains an efficient implementation of the egalitarian mechanism of Bogomolnaia and Moulin.

Comments, suggestions and questions are welcome.

********************************************************************************************************************************************************************

The following files are included:
    - mf.cpp computes a maxmin-fair decomposition
    - bk.cpp computes a list of matchings for each block in the decomposition
    - mf.sh combines both

Example of compiling and running on input "out.flickr-groupmemberships":

    ulimit -s unlimited 
    g++ -O3 -DNDEBUG mf.cpp -o mf && ./mf < sample_graph 

By default mp.cpp outputs a summary of the fair decomposition. To obtain the full decomposition it needs FULL_DECOMP flag needs to be enabled at compile time: 
    g++ -O3 -DNDEBUG -DFULL_OUTPUT mf.cpp -o mf && ./mf < sample_graph | grep '#'  | cut -d' ' -f2- > decomp_file

    g++ -O3 -DNDEBUG bk.cpp -o bk && ./bk decomp_file | grep '#' | cut -d' ' -f2- > matchings_list

An script to perform both operations and include the results and some summary information in an output directory is provided in mf.sh:    
    ./mf.sh sample_graph result_dir/
Note that if a decomposition file exists in result_dir/ (from a previous run of the script), it will not be recomputed. To recompute from scratch, you must erase the contents of
result_dir or use another directory. The script also compiles mf.cpp and bk.cpp if no binary exists in the current directory.

Data is read from standard input. The format for bipartite graphs is a list of pairs of positive integers (i, j) denoting an edge between the ith vertex from the
left side and the jth vertex from the right side. Each input line contains such a pair, and each pair is separated by a space character. An example is provided in
the sample_graph file.

Certain parts of the the maximum-flow code used in this software (specifically BucketList.cpp) have been based on Andrew Goldberg's high-performance implementation 
of the highest-label push-relabel method (HIPR), available from http://www.avglab.com/andrew/soft.html.

*******************************************************************************************************************************************************************

The program is available on "as is" basis. It is not guaranteed to be free of bugs, and we assume no responsibility for any potential problems. 

This software comes with NO WARRANTY, expressed or implied. By way of example, but not limitation, we make no representations of warranties of merchantability or
fitness for any particular purpose or that the use of the software components or documentation will not infringe any patents, copyrights, trademarks, or other
rights.
