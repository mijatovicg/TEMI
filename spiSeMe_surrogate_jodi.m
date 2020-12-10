function ieiSurrogates = spiSeMe_surrogate_jodi(ieiSequence, varargin)
% SPISEME_SURROGATE_JODI Generates surrogates of Inter-Event-Intervals (IEI) sequences.
%
%	ieiSurrogates = SPISEME_SURROGATE_JODI(ieiSequence, ...)
%	Generates surrogate sequences of IEI corresponding to the original
%	sequence stored in the ieiSequence array by means of the
%	JOint DIstribution (JODI) algorithm.
%
%	Options (passed as 'name', value pairs):
%
%	'M'		Number of surrogate sequences to be generated. By
%			default, M = 1. Each of the M surrogate sequences
%			is a column in the returned array, which therefore
%			has size LxM, where L is the original sequence
%			length.
%	
%	'Verbose'	Sets the verbosity of the function. If true (default),
%			all messages are printed on the command line. If
%			false, only critical errors are displayed.
%
%
%	REF:
%	The JOint DIstribution method (JODI) for generation of surrogate
%	event sequences was originally	proposed by L. Ricci et al. in
%	Chaos 29 (2019), 121102, <a href="matlab:web('https://doi.org/10.1063/1.5138250')">doi:10.1063/1.5138250</a>
%
%	This function is part of the SpiSeMe package.


	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'ieiSequence', @isvector);
	addParameter(ip, 'M', 1, @isnumeric);
	addParameter(ip, 'verbose', true, @islogical);
	ip.KeepUnmatched = false;
	parse(ip, ieiSequence, varargin{:});
	if (~iscolumn(ieiSequence))
		ieiSequence = ieiSequence';
	end
	L = length(ieiSequence);
	M = fix(ip.Results.M);
	if (M < 1)
		error('Function argument "M" must be a positive integer.')
	end
	beVerbose = ip.Results.verbose;

	% --- Provide some information
	if (beVerbose)
% 		fprintf('\n### Starting JODI routine ###\n');
	end

	% --- Extract sequence of ranks out of IEI sequence
	[sortedOriginalIEI, orderOriginalIEI] = sort(ieiSequence);
	ranks = (1:L)';
	sequenceOriginalRanks(orderOriginalIEI,1) = ranks;

	% --- Build histograms according to the Freedman-Diaconis rule
	[h1counts, ~] = histcounts(sequenceOriginalRanks, 'BinMethod', 'fd');
	nBins = length(h1counts);
	[h2counts, h2Xedges, h2Yedges] = histcounts2(sequenceOriginalRanks(1:L-1), sequenceOriginalRanks(2:L), [nBins nBins]);
	if (nBins < 2)
		warning('Histogram has only one bin!');
	end

	% --- Initialize output array (each of the M generated surrogates is a column)
	ieiSurrogates = [];

	% --- Generate surrogates: iterate M times
	for iter=1:M

		if(beVerbose)
% 			fprintf('# Surrogate number %d out of %d.\n', iter, M);
		end

		% --- Initialize surrogate ranks generation
		sequenceSurrogateRanks = zeros(L, 1);
		rng('shuffle', 'twister');

		% --- Generate the first pair r_1, r_2 according to the sample joint distribution
		countUpTo = randi(L);
		nCount = 0;
		stopSearching = false;
		for i = 1:nBins
			for j = 1:nBins
				nCount = nCount + h2counts(i, j);
				if (nCount >= countUpTo)
					stopSearching = true;
					break;
				end
			end
			if (stopSearching)
				break;
			end
		end
		sequenceSurrogateRanks(1) = h2Xedges(i) + rand(1)*(h2Xedges(i+1) - h2Xedges(i));
		sequenceSurrogateRanks(2) = h2Yedges(j) + rand(1)*(h2Yedges(j+1) - h2Yedges(j));

		% --- Iterate the Markov chain to generate successive intervals
		n = 2;
		while (n < L)
			% --- Compute the conditional distribution P(r_{n+1} | r_n) corresponding to r_n being the last generated element
			last_extracted_bin = discretize(sequenceSurrogateRanks(n), h2Xedges);
			conditionalDistribution = h2counts(last_extracted_bin, :);
			normConditionalDistribution = sum(conditionalDistribution);

			% --- Generate the new element according to the conditional distribution P(r_{n+1} | r_n)
			countUpTo = randi(normConditionalDistribution);
			nCount = 0;
			for j = 1:nBins
				nCount = nCount + conditionalDistribution(j);
				if (nCount >= countUpTo)
					break;
				end
			end
			sequenceSurrogateRanks(n+1) = h2Yedges(j) + rand(1)*(h2Yedges(j+1) - h2Yedges(j));
			n = n + 1;
		end

		% --- Transform ranks back into IEIs
		[~, orderSurrogateRanks] = sort(sequenceSurrogateRanks);
		thisSurrogate(orderSurrogateRanks, 1) = sortedOriginalIEI;

		if (beVerbose)
% 			fprintf('Process ended; %d elements generated.\n\n', n);
		end

		% --- Append sequence to the set of those already generated
		ieiSurrogates = [ieiSurrogates, thisSurrogate];
	end
end
