% Classification GUI and toolbox
% Version 1.0
%
% GUI start commands
%
%		classifier			- Start the classification GUI
%		enter_distributions - Starts the parameter input screen (used by classifier)
%		multialgorithms		- Start the algorithm comparison screen
%       show_algorithms     - Show possible algorithms for classification, clustering and feature selection
%
% Preprocessing methods
%
%		ADDC				- Compute k clusters for the data using the agglomerative clustering method
%		AGHC				- Compute k clusters for the data using the agglomerative hierarchical clustering method
%		BIMSEC				- Compute k clusters for the data using the basic iterative MSE clustering method
%		Competitive_learning - Compute k clusters for the data using a competitive NN
%		Deterministic_annealing - Compute k patterns which typify the data using the Deterministic SA algorithm
%		Deterministic_SA    - Compute k patterns which typify the data using the Deterministic SA algorithm (Another implementation)
%		DSLVQ				- Distinction sensitive linear vector quantization 
%		Fuzzy_k_means		- Compute k means for the data using the fuzzy_k_means algorithm
%		FishersLinearDiscriminant	- Fisher linear discriminant
%		k_means				- Compute k means for the data using the k_means algorithm
%       Kohonen_SOFM        - Reduce data points using a Kohonen self-orgenizing feature map
%		Leader_Follower		- Compute clusters for the data using the basic leader-follower clustering method
%		LVQ1				- Linear vector quantization with one neighbor
%		LVQ3				- Linear vector quantization 3 algorithm
%       min_spanning_tree   - Reduce data points using a minimum spanning tree
%		PCA					- Principal component analysis
%       SOHC                - Compute k clusters for the data using the stepwise optimal hierarchical clustering method
%		Stochastic_SA       - Compute k patterns which typify the data using the Stochastic SA algorithm
%
% Parametric classification algorithms
%
%       Balanced_Winnow     - Balanced Winnow algorithm
%       Bayesian_Model_Comparison - Find a Gaussian model using Bayesian model comparison
%		EM					- Expectation maximization algorithm
%       Gibbs               - The Gibbs algorithm
%       Ho_Kashyap          - The regular and modified Ho-Kashyap algorithm
%       LMS                 - Least-means square algorithm
%		LS					- Least squares algorithm
%       Marginalization     - Classify when a feature is missing using the marginal distribution
%		ML					- Maximum likelihood algorithm
%		ML_diag	    		- Maximum likelihood with diagonal covariance matrices
%       ML_II               - Find a Gaussian model using maximum likelihood model comparison
%       NDDF                - Normal density discriminant function
%		None					- A dummy file 
%		Perceptron			- Single perceptron algorithm
%		Perceptron_Batch	- Batch perceptron algorithm
%		Perceptron_BVI		- Batch variable increment perceptron algorithm
%		Perceptron_FM		- Perceptron which improves according to the example farthest from the margin
%		Perceptron_VIM		- Variable increment perceptron with margin algorithm
%		Pocket				- Pocket algorithm
%		RDA				    - Regularized descriminant analysis (Friedman shrinkage algorithm)
%       Relaxation_BM       - Batch relaxation with margin
%       Relaxation_SSM      - Single-sample relaxation with margin
%       Stumps              - Simple stump classifier
%
% Non-parametric classification algorithms
%
%       Ada_Boost                   - Ada Boost algorithm
%       Backpropagation_Batch       - Neural network trained with the batch backpropagation algorithm
%       Backpropagation_CGD         - Neural network trained with the batch backpropagation algorithm and conjugate gradient descent
%       Backpropagation_Quickprop   - Neural network trained with the quickprop backpropagation algorithm
%       Backpropagation_Recurrent   - A recurrent neural network trained with the batch backpropagation algorithm
%       Backpropagation_SM			   - A recurrent neural network trained with the stochastic backpropagation algorithm with momentum
%       Backpropagation_Stochastic  - Neural network trained with the stochastic backpropagation algorithm
%       C4_5                        - The C4.5 algorithm
%       CART                        - Classification and regression trees
%       CARTfunctions               - Used by CART
%       Cascade_Correlation         - Cascade-correlation type neural network
%       Components_with_DF          - Component classifiers with descriminant functions
%       Components_without_DF       - Component classifiers without descriminant functions
%       Deterministic_Boltzmann     - Deterministic Boltzmann learning
%       Discrete_Bayes              - Bayes classifier for discrete patterns
%       Genetic_algorithm	        - Basic genetic algorithm
%       Genetic_programming         - Genetic programming of a solution
%       ID3                         - Quinlan's ID3 classification tree algorithm
%       Interactive_Learning        - Interactive learning (Learning with queries)
%		Local_Polynomial	        - Local polynomial fitting
%		loglikelihood	            - Used by Local polynomial fitting
%       LocBoost                    - Local boosting
%       LocBoostFunctions           - Used by LocBoost
%       Minimum_Cost                - Classify under a minimum cost strategy with histogram equalization
%       Multivariate_Splines        - Multivariate adaptive regression splines
%		Nearest_Neighbor	        - Nearest neighbor algorithm
%		NearestNeighborEditing      - Nearest neighbor editing algorithm
%       Optimal_Brain_Surgeon       - Train a backprop. Neural net. and prune it using the optimal brain surgeon algorithm
%		Parzen				        - Parzen window algorithm
%       Perceptron_Voted            - Voted perceptron algorithm.
%		PNN					        - Probabilistic neural network
%       Projection_pursuit          - Projection pursuit regression for classification
%		RCE					        - Reduced coulomb energy algorithm
%       RBF_Network                 - Train a radial-basis function neural network
%		Store_Grabbag		        - An improvement on the nearest neighbor algorithm
%       SVM                         - Support vector machines
%
% Feature selection
%
%       combinations                - Return all combinations of indices (Used by sequential and exhaustive feature selection)
%       Exhaustive_Feature_Selection- Exhaustive feature selection
%       Genetic_Culling             - A Culling type genetic algorithm for feature selection
%       HDR                         - Hierarchical dimensionality reduction
%       ICA                         - Independent component analysis
%		infomat						- Generates the mutual information matrix. Used by Information_based_selection
%		Information_based_selection	- Choose the most relevant features using the Koller-Sawami algorithm
%       MDS                         - Multidimensional scaling
%       MultipleDiscriminantAnalysis- Multiple descriminant analysis
%       NLPCA                       - Non-linear PCA
%       PCA                         - Principle component analysis
%       Sequential_Feature_Selection- Sequential (Forward/Backward) feature selection
%
% Error estimation
%
%		calculate_error		        - Calculates the classification error given a decision surface
%		classification_error        - Used by claculate_error
%		classify_paramteric         - Builds a decision region for multi-Gaussian distributions
%
% Error bounds
%
%		Bhattacharyya 
%		Chernoff
%		Discriminability
%
% GUI housekeeping functions
%
%		calculate_region	        - Finds the data scatter region
%		classifier_commands	        - Classifier screen commands
%		click_points				- Graphically enter a distribution
%		enter_distribution_commands	- Used by enter_distributions
%       feature_selection           - The feature selection GUI open when data with more than 2D is loaded
%       feature_selection_commands  - The commands file for the feature selection GUI
%       FindParameters              - A GUI for finding the optimal parameters for a classifier
%       FindParametersFunctions     - The commands file for FindParameters
%       GaussianParameters          - Opens a GUI for displaying the gaussian parameters of a distribution
%		generate_data_set	        - Generate a data set given Gaussian parameters
%		high_histogram              - Generate a histogram for high-dimensional data
%		load_file			        - Load data files
%		make_a_draw			        - Randomly find indices from a data set
%       multialgorithms_commands    - Multialgorithms screen comands
%       plot_process                - Plot partition centers during the algorithm execution
%		plot_scatter		        - Make a scatter plot of a data set
%       Predict_performance         - Predict performance of algorithms from their learning curves
%       process_params              - Read a parameter vector and return it's components
%		read_algorithms	            - Reads an algorithm file into a data structure
%		start_classify		        - Main function used by classifier
%		voronoi_regions		        - Plot Voronoi regions
%
% Data sets (Ending _data means that the file contains patterns, 
%                   _params means that the file contains the distribution parameters)
%
%		chess                       - The parameters for a 4x4 chess board distribution
%		clouds				        - A data set composed of four Gaussians 
%		seperable			        - A linearly seperable data set
%		spiral				        - Two interlocking spirals data set
%		XOR							- XOR distribution
%
%
%____________________________________________________________________________________
%  Elad Yom-Tov (elad@ieee.org) and David Stork
%  Technion - Israel Institute of Technology
%  Haifa, Israel