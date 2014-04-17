import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.HashSet;



public class PathExpressionComparator implements Comparator<Object> {


	HashMap<PairPath,HashSet<List<Integer>>> pp_to_transcript;
	HashMap<List<Integer>,Float> transcript_to_expr;
	HashMap<List<Integer>,HashMap<PairPath,Float>> transcript_to_pp_frag_counts; // split out by pairpath
	HashMap<List<Integer>, Float> transcript_to_sum_frag_counts; // sum of all fractional frag counts
	
	HashMap<PairPath,Integer> pp_to_read_support;
	
	
	public static int MAX_ROUNDS = 100;
	public static int MIN_ROUNDS = 100;
	public static float MIN_DELTA_REQUIRED = 0.01f;
	
	
	
	
	
	public PathExpressionComparator (List<List<Integer>> all_paths, HashMap<List<Integer>,HashMap<PairPath, Integer>> PathReads) {
		
		pp_to_transcript = new HashMap<PairPath,HashSet<List<Integer>>>();
		transcript_to_expr = new HashMap<List<Integer>,Float>();
		transcript_to_pp_frag_counts = new HashMap<List<Integer>,HashMap<PairPath,Float>>();
		pp_to_read_support = new HashMap<PairPath,Integer>();
		transcript_to_sum_frag_counts = new HashMap<List<Integer>, Float>();
		
		BFLY_GLOBALS.debugMes("EM on: " + all_paths.size() + " paths.", 10);
		
		//if (VERBOSE)
		//	System.err.println("-reorganizing data structures for EM computation.");
		
		// initialize
		for (List<Integer> path : PathReads.keySet()) {
			
			
			// only analyze the paths of interest, not all stored in the more comprehensive PathReads data structure.
			if (! all_paths.contains(path))
				continue;
			
			
			transcript_to_pp_frag_counts.put(path, new HashMap<PairPath,Float>());
			transcript_to_sum_frag_counts.put(path, 0f);
			
			
			HashMap<PairPath,Integer> pp_map = PathReads.get(path);
			pp_to_read_support.putAll(pp_map);
			
			for (PairPath pp : pp_map.keySet()) {
				
				//if (! pp.isCompatibleAndContainedBySinglePath(path))
				//	continue;
				
				if (! pp_to_transcript.containsKey(pp)) {
					pp_to_transcript.put(pp, new HashSet<List<Integer>>());
				}
				pp_to_transcript.get(pp).add(path);
				transcript_to_pp_frag_counts.get(path).put(pp, 0f); // just init for now, store fractional assignments later on.
			}
			
			
		}
		
	
		
		run_EM();
		
		
		
	}
	
	
	
	
	private void run_EM () {
		
		init_transcript_expr();
		
		boolean enough_iterations = false;
		float prev_likelihood = calc_likelihood();
	
		//describe_expr_values();
	
		int round = 0;
		
		while (! enough_iterations) {
			
			round++;
			
			long start_time = System.currentTimeMillis();
			
			long start_time_E_step = System.currentTimeMillis();
			calc_expected_read_alignments ();
			long end_time_E_step = System.currentTimeMillis();
			//System.err.println("Time for E-step: " + (end_time_E_step - start_time_E_step));
			
			long start_time_M_step = System.currentTimeMillis();
			recompute_expression_values();
			long end_time_M_step = System.currentTimeMillis();
			//System.err.println("Time for M-step: " + (end_time_M_step - start_time_M_step));
			
			// describe_expr_values();
			
			long start_time_likelihood = System.currentTimeMillis();
			float likelihood = calc_likelihood();
			long end_time_likelihood = System.currentTimeMillis();
			//System.err.println("Time for likelihood calc: " + (end_time_likelihood - start_time_likelihood));
			
			float delta = likelihood - prev_likelihood;
			
			long end_time = System.currentTimeMillis();
			
			long seconds_for_computation = (end_time - start_time) / 1000l;
			
			BFLY_GLOBALS.debugMes("EM round[" + round + "] Prev: " + prev_likelihood + ", curr: " + likelihood 
						+  ", delta: " + delta + " [seconds: " + seconds_for_computation + "]\n", 10);
		
			prev_likelihood = likelihood;
			
			
			
			if ( Math.abs(delta) < MIN_DELTA_REQUIRED || round >= MAX_ROUNDS)
				break; 
		}
		
	
		
	}
	
	
	private void init_transcript_expr() {
		
		//if (VERBOSE)
		//	System.err.println("// initializing EM (round 0)");
		
		
		init_transcript_to_sum_frags();
		
		
		// loop through pair paths, assign  read_support/num_transcript to each transcript
		for (PairPath pp : pp_to_transcript.keySet()) {
			
			HashSet<List<Integer>> transcript_list = pp_to_transcript.get(pp);
			int num_transcripts_with_pp = transcript_list.size();
			int read_support = pp_to_read_support.get(pp);
			
			// initially, divide multiply mappers among the transcripts
			// including the number of reads in that pp
			float fractional_support = (float)read_support/num_transcripts_with_pp;
			
			for (List<Integer> transcript : transcript_list) {
				
				// add fractional support to the expression value of that transcript.
				
				transcript_to_sum_frag_counts.put(transcript, transcript_to_sum_frag_counts.get(transcript) + fractional_support);
				
				
				// store fractional support for this transcript here.
				this.transcript_to_pp_frag_counts.get(transcript).put(pp, fractional_support);  
				
			}
		}
		
		
		// set initial expression values.
		float sum_expr_values = 0;
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			float frags = transcript_to_sum_frag_counts.get(transcript);
			
	
			
			float expr = compute_expr(transcript, frags); 
			
			this.set_expr(transcript, expr);
			
			sum_expr_values += expr;
		}
		
		// convert to relative expression value:
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			float rel_expr_val = get_expr(transcript) / sum_expr_values;
			set_expr(transcript, rel_expr_val);
		}
		
		
		
	}
	
	
	private void init_transcript_to_sum_frags() {
		for (List<Integer> path : transcript_to_sum_frag_counts.keySet()) {
			transcript_to_sum_frag_counts.put(path, 0f);
		}
		
	}




	public String describe_frag_assignment (List<Integer> transcript) {
		
		String text = "";
		
		float sum_frags_assigned = 0;
		for (PairPath pp : transcript_to_pp_frag_counts.get(transcript).keySet()) {
			
			float orig_pp_read_support = this.pp_to_read_support.get(pp);

			float frags_assigned = transcript_to_pp_frag_counts.get(transcript).get(pp);
			
			float pct_pp_reads_assigned = frags_assigned / orig_pp_read_support * 100;
			
			text += "FRAGS_ASSIGNED: [" + frags_assigned + ", " + pct_pp_reads_assigned + "%] for read: " + pp + "\n";
			sum_frags_assigned += frags_assigned;
		}
		
		text = "EM_RESULT: len=" + transcript.size() + " expr=" + get_expr(transcript) + " sum_frags=" + sum_frags_assigned + "\n" + text;
		
		return(text);
		
	}
	
	
	public float get_expr (List<Integer> path) {
		return(transcript_to_expr.get(path));
	}
	
	
	private void describe_expr_values () {
		
		System.err.println("####  ACCOUNTING FOR EXPR VALUES");
		
		float sum_expr = 0;
		
		for (List<Integer> transcript : transcript_to_expr.keySet()) {
			
			float expr = get_expr(transcript);
			sum_expr += expr;
			
			System.err.println("expr: " + expr + " for path: " + transcript);
			
		}
		
		System.err.println("\n** Sum expr: " + sum_expr + "\n\n");	
		
	}
	
	
	
	private Float calc_likelihood () {
		
		//if (VERBOSE)
		//	System.err.println("-computing likelihood");
		
		// compute probability of all observed read alignments, given isoform relative expression values.
		
		//  formula: prod_for_all_reads (sum_all_isoform_mappings ( isoform_relative_expression *  1/length(isoform)) )
		
		float log_likelihood = 0;
		
		for (PairPath pp : pp_to_transcript.keySet()) {
			HashSet<List<Integer>> transcripts = pp_to_transcript.get(pp);
			
			float read_prob = 0;
			for (List<Integer> transcript : transcripts) {
				//int num_pairpaths_in_transcript = transcript_to_pp_frag_counts.get(transcript).size();
				
				float expr = get_expr(transcript);
				
				//read_prob += expr/num_pairpaths_in_transcript;
				
				read_prob += expr;
				
			}
			
			log_likelihood += Math.log(read_prob);
				
		}
	
		
		return(log_likelihood);
		
	}
	
	
	//--------------
	// The E-step
	//--------------
	
	private void calc_expected_read_alignments () {
		
		//if (VERBOSE)
		//	System.err.println("-calc expected read alignments (E-step)");
		
		// fractional read assignment
		
		//  Expected value for read alignment =   prob of read alignment to transcript_i /  sum of prob of read aligned to all other transcripts
		
		//                                    =   (expression(transcript_i) / length(transcript_i) ) / sum_for_all_multiply_mapped_read (expr(trans_o)/length(trans_o)))
		
		
		
		this.init_transcript_to_sum_frags();
		
		//float pp_counter = 0;
		for (PairPath pp : pp_to_transcript.keySet()) {
			HashSet<List<Integer>> transcripts = pp_to_transcript.get(pp);
			
			//pp_counter++;
			//System.err.println("\r[" + pp_counter + "]  ");
			
			
			int count_reads_in_pp = pp_to_read_support.get(pp);
			
			float sum_read_probs = 0;
			for (List<Integer> transcript : transcripts) {
			
				float expr = get_expr(transcript);
				sum_read_probs += expr; 
				
			}
			
			for (List<Integer> transcript : transcripts) {
				
				float expr = get_expr(transcript);
				float this_read_prob = expr;
				
			
				
				float expected_fractional_assignment = this_read_prob / sum_read_probs;
				
				float fractional_frags_assigned = count_reads_in_pp * expected_fractional_assignment;
				
			
				transcript_to_sum_frag_counts.put(transcript, transcript_to_sum_frag_counts.get(transcript) + fractional_frags_assigned);
				
				// store fractional support for this transcript here.
				this.transcript_to_pp_frag_counts.get(transcript).put(pp, fractional_frags_assigned); 
			}
			
		}
		

		
	}
	
	//------------
	// The M-step
	//------------
	
	private void recompute_expression_values () {
	
		//if (VERBOSE)
		//	System.err.println("-recomputing fractional expression values (M-step)");
		
		float sum_expr_vals = 0;
		
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			
			float frag_counts = transcript_to_sum_frag_counts.get(transcript);
			//int num_pairpaths_in_transcript = transcript_to_pp.get(transcript).size();
			
			float expr = compute_expr(transcript, frag_counts); // / num_pairpaths_in_transcript;
		
			// temporary replace w/ expr value
			set_expr(transcript, expr);
			sum_expr_vals += expr;
			
		}
		
		// update w/ fractional expression value
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			
			float expr = get_expr(transcript);
			
			float fractional_expression = expr/sum_expr_vals;
			
			set_expr(transcript, fractional_expression);
				
			
		}
		
	}
	

	private float compute_expr(List<Integer> transcript, float frag_counts) {
		
		return(frag_counts);
		
	}
	
	
	private void set_expr(List<Integer> transcript, float fractional_expression) {
		
		this.transcript_to_expr.put(transcript, fractional_expression);
		
	}




	@Override
	public int compare(Object o1, Object o2) {
		
		List<Integer> path1 = (List<Integer>) o1;
		List<Integer> path2 = (List<Integer>) o2;
		
		Float expr_path1 = transcript_to_expr.get(path1);
		Float expr_path2 = transcript_to_expr.get(path2);
		
		if (expr_path1 < expr_path2) {
			return(-1);
		}
		else if (expr_path1 > expr_path2) {
			return(1);
		}
		else {
			return(0);
		}
	
	}
	
	
		
}
	
	
	
