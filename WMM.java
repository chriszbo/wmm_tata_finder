import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.List;
import java.util.ArrayList;

public class WMM {
	private static final int UNMAPPED = 4;
	private static final int LOWQUALITY = 512;
	private static final Pattern POLYA = Pattern.compile("([ARWMDHVN]+)$");

	/*   1 2 3 4 5 6
	 * A
	 * C
	 * G
	 * T
	 */
	private static final int A = 0;
	private static final int C = 1;
	private static final int G = 2;
	private static final int T = 3;

    // Format is whole number int percentages
	private static final int W0[][] = {{100, 100, 0, 100, 100, 100}, 
	                                  {0, 0, 0, 0, 0, 0}, 
	                                  {0, 0, 0, 0, 0, 0}, 
	                                  {0, 0, 100, 0, 0, 0}};

	private static final int W1[][] = {{85, 85, 5, 85, 85, 85}, 
	                                  {5, 5, 5, 5, 5, 5}, 
	                                  {5, 5, 5, 5, 5, 5},
	                                  {5, 5, 85, 5, 5, 5}};

	private static final double BACKGROUND_A = 0.293;
	private static final double BACKGROUND_G = 0.207;
	private static final double BACKGROUND_C = 0.200;
	private static final double BACKGROUND_T = 0.300;

    // Current constants result in 1088 filtered reads.
	public static final int MAXALIGNMENTSCORE = -3;
	public static final int MINMISMATCHES = 5;
	public static final int MINPOLYALENGTH = 10;
	public static final int MAXPOLYALENGTH = 55;

    // TODO: Return a new weight matrix formed by performing meme algo on weightMatrix.
	public static int[][] meme(int[][] weightMatrix) {
		return null;
	}

    // TODO: Return an array of normalizing terms.
    // Get denominators for each column to normalize the weights. Each normalizing term will be
	// the sum of all the weights in each column.
	public static int[] getNormalizingTerms(int[][] weightMatrix) {
		return null;
	}

	public static void main(String[] args) {
		// use backslash on windows, forward slash on linux
	 	String path = System.getProperty("user.dir") + "/all.sam";
	 	List<FilteredRead> filteredReads = null;
	 	try {
	 		filteredReads = WMM.filterReads(path);
	 	} catch (IOException e) {
	 		e.printStackTrace();
	 	}
	 	int defaultNormalizingTerm[] = {100, 100, 100, 100, 100, 100};
	 	ScoringReport s0 = scoreSequences(filteredReads, W0, defaultNormalizingTerm);
	 	ScoringReport s1 = scoreSequences(filteredReads, W1, defaultNormalizingTerm);

	 	// calculate W2 with meme algo
	 	int W2[][] = meme(W1);
	 	//
	 	int[] normalizingTerm = getNormalizingTerms(W2);

	 	ScoringReport s2 = scoreSequences(filteredReads, W2, normalizingTerm);
	 	System.out.println("Model: WMM0");
	 	System.out.println("Number of Candidates with Positive LLR: " + s0.getCandidatesPositiveLLR());
	 	System.out.println("Distance between cleavate site and left end of best scoring hit: " + s0.getDistanceTATAToCleavage());
	 	System.out.println("Relative Entropy: " + s0.getRelativeEntropy());
	 	System.out.println("-------------------------------------------");
	 	System.out.println("Model: WMM1");
	 	System.out.println("Number of Candidates with Positive LLR: " + s1.getCandidatesPositiveLLR());
	 	System.out.println("Distance between cleavate site and left end of best scoring hit: " + s1.getDistanceTATAToCleavage());
	 	System.out.println("Relative Entropy: " + s1.getRelativeEntropy());
	 	System.out.println("-------------------------------------------");
	 	System.out.println("Model: WMM2");
	 	System.out.println("Number of Candidates with Positive LLR: " + s2.getCandidatesPositiveLLR());
	 	System.out.println("Distance between cleavate site and left end of best scoring hit: " + s2.getDistanceTATAToCleavage());
	 	System.out.println("Relative Entropy: " + s2.getRelativeEntropy());
	 }

	 public static ScoringReport scoreSequences(List<FilteredRead> filteredReads, 
	 	                                        int[][] weightMatrix, int[] normalizingTerm) {
	 	int candidatesPositiveLLR = 0;
	 	int sumDistance = 0;
	 	for (FilteredRead read : filteredReads) {
	 		String sequence = read.getSequence();
	 		int cleavageIndex = read.getCleavageIndex();
	 		if (cleavageIndex < 6) {
	 			continue;
	 		}
	 		double maxLLR = Double.NEGATIVE_INFINITY;
	 		int distance = 0;
	 		for (int i = 0; i < cleavageIndex - 5; i++) {
	 			String frame = sequence.substring(i, i + 6);
	 			double LLR = scoreFrame(frame, weightMatrix, normalizingTerm);
	 			if (LLR > maxLLR) {
	 				maxLLR = LLR;
	 				distance = read.getCleavageIndex() - i;
	 			}
	 		}
	 		if (maxLLR > 0) {
	 			candidatesPositiveLLR++;
	 			sumDistance += distance;
	 		}
	 	}
	 	double avgLLR = (double) sumDistance / (double) candidatesPositiveLLR;
	 	// calculate relative entropy
	 	double relEntropy = relativeEntropy(weightMatrix, normalizingTerm);
	 	return new ScoringReport(candidatesPositiveLLR, avgLLR, relEntropy);
	 }

	 public static double relativeEntropy(int[][] weightMatrix, int[] normalizingTerm) {
	 	double wholeSum = 0;
	 	for (int i = 0; i < weightMatrix[0].length; i++) {
	 		double colSum = 0;
	 		for (int j = 0; j < weightMatrix.length; j++) {
	 			if (weightMatrix[j][i] == 0) {
	 				continue;
	 			}
	 			double ratio;
	 			if (j == A) {
	 				ratio = (double) weightMatrix[j][i] / BACKGROUND_A;
	 			} else if (j == C) {
	 				ratio = (double) weightMatrix[j][i] / BACKGROUND_C;
	 			} else if (j == G) {
	 				ratio = (double) weightMatrix[j][i] / BACKGROUND_G;
	 			} else {
	 				ratio = (double) weightMatrix[j][i] / BACKGROUND_T;
	 			}
	 			ratio /= normalizingTerm[i];
	 			colSum += (((double) weightMatrix[j][i] / normalizingTerm[i]) * (Math.log(ratio) / Math.log(2)));
	 		}
	 		wholeSum += colSum;
	 	}
	 	return wholeSum;
	 }

	 public static double scoreFrame(String frame, int[][] weightMatrix, int[] normalizingTerm) {
	 	double sum = 0;
	 	for (int i = 0; i < 6; i++) {
	 		char c = frame.charAt(i);
	 		int weight;
	 		double ratio;
	 		if (c == 'A') {
	 			weight = weightMatrix[A][i];
	 			if (weight == 0) {
	 				return Double.NEGATIVE_INFINITY;
	 			} else {
	 				ratio = (double) weight / BACKGROUND_A;
	 			}
	 		} else if (c == 'C') {
	 			weight = weightMatrix[C][i];
	 			if (weight == 0) {
	 				return Double.NEGATIVE_INFINITY;
	 			} else {
	 				ratio = (double) weight / BACKGROUND_C;
	 			}
	 		} else if (c == 'G') {
	 			weight = weightMatrix[G][i];
	 			if (weight == 0) {
	 				return Double.NEGATIVE_INFINITY;
	 			} else {
	 				ratio = (double) weight / BACKGROUND_G;
	 			}
	 		} else {
	 			weight = weightMatrix[T][i];
	 			if (weight == 0) {
	 				return Double.NEGATIVE_INFINITY;
	 			} else {
	 				ratio = (double) weight / BACKGROUND_T;
	 			}
	 		} 
	 		sum += Math.log(ratio / normalizingTerm[i]) / Math.log(2);
	 	}
	 	return sum;
	}

	public static List<FilteredRead> filterReads(String fileName) throws IOException {
		FileInputStream fstream = new FileInputStream(fileName);
		BufferedReader reader = new BufferedReader(new InputStreamReader(fstream));
		List<FilteredRead> sequences = new ArrayList<FilteredRead>();
		// consume lines until start of sequence data
		while (!reader.readLine().startsWith("@PG")) {}
		int count = 1;
	    int skip1 = 0;
	    int skip2 = 0;
	    int skip3 = 0;
	    int skip4 = 0;
	    int skip5 = 0;
	    int total = 0;
		while (true) {
			String line = reader.readLine();
			if (line == null) {
				break;
			}
			total++;
			String[] fields = line.trim().split("\\s+");
			int flag = Integer.parseInt(fields[1]);
			String cigar = fields[5].trim();
			String readID = fields[0];
			String location = fields[2].trim() + " " + fields[3].trim();
			String sequence = fields[9].trim();
			if (flag == UNMAPPED || flag == LOWQUALITY) {
				skip1++;
				continue;
			}
			int alignmentScore = Integer.parseInt(fields[11].trim().substring(5));
			if (alignmentScore > MAXALIGNMENTSCORE) {
				skip2++;
				continue;
			}
			int mismatches = Integer.parseInt(fields[16].trim().substring(5));
			if (mismatches < MINMISMATCHES) {
				skip3++;
				continue;
			}
			Matcher matcher = POLYA.matcher(sequence);
			if (matcher.find()) {
				if (matcher.group().length() < MINPOLYALENGTH || matcher.group().length() > MAXPOLYALENGTH) {
					// poly A tail length too small or too long
					skip4++;
					continue;
				}
			} else {
				// Poly A tail not found
				skip5++;
				continue;
			}
			count++;
			// replace non ATGC in sequence with T?
			sequence = sequence.replaceAll("[^ATGCatgc]", "T").toUpperCase();
			sequences.add(new FilteredRead(readID, location, sequence, matcher.start()));
		}
		System.out.println("Total Reads Examined: " + total);
		System.out.println("Candidates Selected: " + count);
		System.out.println("skip1 omits: " + skip1);
		System.out.println("skip2 omits: " + skip2);
		System.out.println("skip3 omits: " + skip3);
		System.out.println("skip4 omits: " + skip4);
		System.out.println("skip5 omits: " + skip5);
		return sequences;
	}

	public static class ScoringReport {
		private int candidatesPositiveLLR;
		private double distanceTATAToCleavage;
		private double relativeEntropy;
		public ScoringReport(int candidatesPositiveLLR, 
			                 double distanceTATAToCleavage, 
			                 double relativeEntropy) {
			this.candidatesPositiveLLR = candidatesPositiveLLR;
			this.distanceTATAToCleavage = distanceTATAToCleavage;
			this.relativeEntropy = relativeEntropy;
		}
		public int getCandidatesPositiveLLR() {
			return candidatesPositiveLLR;
		}
		public double getDistanceTATAToCleavage() {
			return distanceTATAToCleavage;
		}
		public double getRelativeEntropy() {
			return relativeEntropy;
		}
	}

	public static class FilteredRead {
		private String readID;
		private String location;
		private String sequence;
		private int cleavageIndex;
		public FilteredRead(String readID, String location, String sequence, int cleavageIndex) {
			this.readID = readID;
			this.location = location;
			this.sequence = sequence;
			this.cleavageIndex = cleavageIndex;
		}
		public String getReadID() {
			return readID;
		}
		public String getLocation() {
			return location;
		}
		public String getSequence() {
			return sequence;
		}
		public int getCleavageIndex() {
			return cleavageIndex;
		}
		@Override
	    public int hashCode() {
	        return readID.hashCode();
	    }
	    @Override
		public boolean equals(Object obj) {
		    if (obj == null || 
		    		!FilteredRead.class.isAssignableFrom(obj.getClass())) {
		        return false;
		    }
			final FilteredRead other = (FilteredRead) obj;
		    if (readID.equals(other.getReadID())) {
		    	return true;
		    }
		    return false;
		}
	}
}