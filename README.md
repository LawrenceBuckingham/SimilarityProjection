<h1>Similarity Projection and Related Programs</h1>
<h2>Contents</h2>
<ul>
	<li><a href="#Desription">Description</a></li>
	<li><a href="#Publications">Publications</a></li>
	<li><a href="#Requirements">Requirements</a></li>
	<li>Get started with Similarity Projection and Compressed <em>k</em>-mer Vocabularies</li>
	<li>Utility programs</li>
	<li>Experimental GUI programs</li>
</ul>
<a name="Description"><h2>Description</h2></a>
<p>This repository contains the current release(s) of the software developed in the course of my PhD candidature, with subsequent revisions.
</p>
<p>
    The main outcome of the research is the Compressed <em>k</em>-mer Vocabulary metric embedding procedure for sequential symbolic data. This method projects (biological) sequences into a metric space by using a technique similar to Vector Quantisation to generate a sparse high-dimensional binary feature vector which represents each sequence. A standard metric in the binary embedding space is then used to estimate the distance between the original sequences.</p>

<p>Work to date primarily focusses on biological sequence data, including protein and RNA, but the approach should be equally applicable to any sequentially structured symbolic data.</p>
<p>The programs included in this repository implement algorithms for alignment-free sequence comparison which have been demonstrated to achieve accuracy superior to that of <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi">NCBI BLAST</a>, with execution speed comparable to that of <a href="http://drive5.com">USEARCH</a>.</p>

<a name="Publications"><h2>Publications</h2></a>
<p>The Similarity Projection concept, algorithm, and preliminary results were published at the 2017 IEEE eScience Conference:</p>
<ul>
    <li>Buckingham, L., Chappell, T., Hogan, J. M., & Geva, S. (2017, 24-27 Oct. 2017). Similarity Projection: A Geometric Measure for Comparison of Biological Sequences. Paper presented at the 2017 IEEE 13th International Conference on e-Science (e-Science).</li>
</ul>
<p>
    The Compressed Kmer Vocabulary metric embedding method was introduced at the 2018 Australasian Document Computing Symposium:
</p>
<ul>
    <li>Buckingham, L., Geva, S., & Hogan, J. M. (2018). Protein database search using compressed k-mer vocabularies. Paper presented at the Proceedings of the 23rd Australasian Document Computing Symposium, Dunedin, New Zealand. https://doi.org/10.1145/3291992.3291997</li>
</ul>
<p>An extended empirical study is under preparation at present; results and associated software will be posted here in due course. This work was severely impacted by workplace disruption associated with COVID-19.</p>

<a name="Requirements">
	<h2>Requirements</h2>
</a>

<p>The programs are written in C++ and compiled with GCC-10, targeting the current 64-bit version of the Cygwin Posix-like environment on Microsoft Windows, and the current Ubuntu flavour of Linux, although I have found that statically linked Ubuntu executables also work fine on Red Hat Enterprise Linux v7.</p>
<p>Source&nbsp;code is mainly edited and compiled with Microsoft Visual Studio 2019 (either the Enterprise or Community Edition should suffice), however all build tasks can be executed from the shell. I have also found Visual Studio Code to be a useful cross-platform alternative with nice integration of GDB, which can be helpful at times for debugging.</p>
<p>Some of these programs use the <a href="https://www.fltk.org/software.php">FLTK library</a> to implement an experimental graphical user interface.</p>

<ul>
    <li>Under Cygwin, use the latest version of FLTK available via the standard Cygwin install/update program.
</li>
    <li>Under Ubuntu, use the command <code>apt install libfltk1.3-dev</code> to install the library.</li>
    <li>To deploy a binary built under Ubuntu to a Linux environment where a suitable compiler is not available, the FLTK object library can be installed on the target machine by:
       
        <ul>
            <li>Create a directory on the target machine. For concreteness, assume the name of the directory is <code>~/libfltk</code>.</li>
            <li>Copy the following files from Ubuntu to the directory on the target machine:<ul>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk.a</code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_gl.so.1.3</code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk.so </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_gl.a </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_cairo.so.1.3 </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_forms.a </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_images.so.1.3 </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_cairo.so </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_images.so </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_images.a </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk.so.1.3 </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_forms.so.1.3 </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_cairo.a </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_gl.so </code></li>
                <li><code>/usr/lib/x86_64-linux-gnu/libfltk_forms.so </code></li>
            </ul>
            </li>
            <li>Add the name of this directory to the dynamic linking search path on the target machine. To do this, edit <code>~/.bashrc</code>, adding the instruction:<br /><code>export LD_LIBRARY_PATH=~/libfltk</code>
            </li>
        </ul>
    </li>
</ul>
<a name="#Get_started"><h2>Get started with Similarity Projection and Compressed <em>k</em>-mer Vocabularies</h2></a>
<p>
	Download (and if necessary extract) the source tree, placing it in a convenient location in the local file system. The directory should resemble the following </p>
<p>
	&nbsp;</p>

<h2>Utility programs</h2>
<p>&nbsp;</p>
<h2>Experimental GUI programs</h2>
<h2>Work in progress...</h2>

I'm still working on this document!
Content will be uploaded in the near future.

