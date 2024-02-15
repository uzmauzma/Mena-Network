<h1><b>Molecular Ecological Network Analysis (MENA)</b></h1>
<h2>Introduction</h2>
<p>This repository contains code related to the investigation of Eukaryote-Prokaryote interactions through Molecular Ecological Network Analysis (MENA) and network analysis techniques. MENA, available at http://129.15.40.240/mena/, is a computational framework commonly used in microbial ecology to examine complex relationships and interactions between taxa of microbial profiles.</p>

<h2>Dataset</h2>
<p>The MENA algorithm was applied to a combined dataset representing the relative abundance of species across different depths (2 cm, 6 cm, 8 cm, 12 cm, and 23 cm) and time points (5 weeks, 9 weeks, 12 weeks, and 23 weeks), encompassing 75% of the samples. The dataset comprised both Eukaryote and Prokaryote data. </p>

<h2>Code Usage Instructions</h2>
<h3>Preparing Your Data for MENA (1_Mena_Readyfile.R)</h3>
<p>This script is used to prepare your data for analysis with MENA. Follow these steps to use it:</p>
<ol>
  <li>The script will preprocess your data and generate files necessary for MENA analysis.</li>
  <li>If the generated dataset containing the top 75% most abundant taxa has fewer than 50 taxa, MENA will not be applied, and the threshold will be adjusted accordingly.</li>
  <li>Adjustments to the minimum prevalence are based on the dataset size to ensure effective analysis with MENA.</li>
</ol>

<h3>Applying MENA to Your Dataset</h3>
<p>Once your data is prepared, you can apply MENA using the following steps:</p>
<h4>1. Upload Datasets:</h4>
<ul>
  <li>Upload the prepared datasets into the MENA analysis environment or software.</li>
</ul>

<h4>2. Construct Network:</h4>
<ul>
  <li>Set the majority parameter to 1 and apply the UTM to detect the threshold.</li>
  <li>Download the two generated files: Raw similarity matrix and condones out table (for the names of the taxa order in the correlation matrix).</li>
  <li>Ensure that other parameters remain the same.</li>
</ul>
<p>Upon completion of this step, two files will be generated: one showing the 'Pearson_Correlation.txt' and the other 'MV_Estimated.txt'. These two files will be used as input for the script "2_Generate_UTM.R" to generate the file 'Upper_triangular_mat.txt', which, in turn, will be used in script "3_Network.R" to generate the "Network.csv" file.</p>

<h4>3. Analyze Network:</h4>
<ul>
  <li>Analyze the network for global network properties.</li>
  <li>Analyze individual node certainty .</li>
  <li>Perform network module separation and modularity calculation (e.g., using fast_greedy modularity algorithm).</li>
</ul>
<p>Upon completion of this step, two files will be generated: one showing the "centrality.txt" and the other "fast_greedy_modularity.txt". These two files will be used as input for the script "4_Node_class.R" to generate the file "Node_attribute.txt". </p>

<h4>4. Visualize with Cytoscape:</h4>
<p>Import the two below files into Cytoscape software.</p>
<ul>
  <li><b>1. "Node_attribute.txt" </b></li>
  <li><b>2. "Network.csv" </b></li>
</ul>
<li>Use Cytoscape's visualization tools to create visual representations of the network.</li>
