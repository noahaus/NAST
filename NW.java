import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class NW {

	public static void main(String[] args) throws IOException, InterruptedException {
		long startTime = System.nanoTime();
		
		/*
		 * READ IN DATA
		 */
		
		FileReader seq1 = new FileReader("C:\\Users\\Owner\\Documents\\workspace\\Needleman_Wunsch\\src\\FastaSampleFiles\\HIV1a.fasta");
		FileReader seq2 = new FileReader("C:\\Users\\Owner\\Documents\\workspace\\Needleman_Wunsch\\src\\FastaSampleFiles\\HIV1b.fasta");
        BufferedReader seq1_buff = new BufferedReader(seq1);
        BufferedReader seq2_buff = new BufferedReader(seq2);
   
        /*
		 * USE THREADS TO READ IN DATA CONCURRENTLY. WAIT TILL LOADED
		 */
        
        InputThread it1 = new InputThread(seq1_buff);
        InputThread it2 = new InputThread(seq2_buff);
        it1.start();
        it2.start();
        it1.join();
        it2.join();
        
        /*
		 * PERFORM NEEDLEMAN-WUNSCH GLOBAL ALGORITHM. DETAILS IN METHOD
		 */
        
        String seq1_final = it1.getFin();
        String seq2_final = it2.getFin();
        long read_data_time = System.nanoTime();
        System.out.println("Time to read in data: "+(read_data_time - startTime)/1000000000.0 + " seconds");
        
		Aligner align = new Aligner(seq1_final,seq2_final);
		long endTime = System.nanoTime();
		System.out.println("Time to do entire algorithm: "+(endTime - startTime)/1000000000.0 + " seconds");
	}

}

class Aligner{
	//Our DNA/RNA/protein sequence
	String seq1="", seq2="";
	//Scores
	int match=1, mismatch=-1, indel=-1;
	//Best alignment string
	String out_seq1="", out_seq2="",buffer="";
	//How we find the best alignment 
	int[][] score_matrix;
	String[][] trace_matrix;
	
	//Constructor, User is okay with defaults
	public Aligner(String seq1, String seq2) {
		this.seq1 =seq1; this.seq2 = seq2;
		//Align now
		align(this.seq1,this.seq2);
	}
	
	//Constructor, User has different scoring system
	public Aligner(String seq1, String seq2, int match, int mismatch, int indel) {
		this.seq1 = seq1; this.seq2 = seq2;
		this.match = match; this.mismatch = mismatch; this.indel = indel;
		
		//Align now
		align(this.seq1,this.seq2);
	}
	
	//
	private void align(String seq1, String seq2) {
		//Initialize the score & trace matrix
		//quicker just to initialize them both in one loop
		char[] seq1_char = seq1.toCharArray();  char[] seq2_char = seq2.toCharArray();
		long align_start = System.nanoTime();
		score_matrix = new int[seq1.length()+1][seq2.length()+1];
		trace_matrix = new String[seq1.length()+1][seq2.length()+1];
		
		int diag = 0;
		int top  = 0;
		int left = 0;
		int temp = 0;
		int max = 0;
		int diag_score = 0;
		int top_score = 0;
		int left_score = 0;
		int[] score_list = new int[3];
		
		int is1,is2 = 0; 
		
		//How to do this faster?
		//if-else versus switch-case
		//switch case wins!
		
		for(int i = 0;i<score_matrix.length;i++) {
			for(int j = 0;j<score_matrix[0].length;j++) {
				
				is1 = i == 0 ? 1:0;
				is2 = j == 0 ? 1:0;
				
				switch(is1) {
					case 1:
						switch(is2) {
						case 1:
							
							trace_matrix[i][j] = "done";
							score_matrix[i][j] = 0;	
							break;
							
						case 0:
							
							score_matrix[i][j] = -j;
							trace_matrix[i][j] = "left";
							break;
						}
					break;
					case 0:
						switch(is2) {
						case 1:
							
							score_matrix[i][j] = -i;
							trace_matrix[i][j] = "up";
							break;
							
						case 0:
							
							diag = score_matrix[i-1][j-1];
							top = score_matrix[i][j-1];	
							left = score_matrix[i-1][j];
							
						
							if(seq1_char[i-1] == seq2_char[j-1]) {
								diag_score = diag + match;
							}
							else {
								diag_score = diag + mismatch;
							}
							
							top_score = top + indel;
							left_score = left + indel;
							
							score_list[0] = diag_score;
							score_list[1] = top_score;
							score_list[2] = left_score;
							
							temp = find_max(score_list);
							score_matrix[i][j] = temp;
							
							 max = temp;
							
							if(max == diag_score)
								trace_matrix[i][j] = "diag";
							else if(max == top_score)
								trace_matrix[i][j] = "left";
							else if(max == left_score)
								trace_matrix[i][j] = "up";
							
							break;
						}
					break;
				}
				
				
			}
			
		}
		
		long align_fill_matrix = System.nanoTime();
		
		System.out.println("Time to fill in the matrices: "+(align_fill_matrix - align_start)/1000000000.0 + " seconds");
		
		/*
		 * MATRICIES FILLED, NOW PRINT THE ALIGNMENT
		 */
				
		this.best_alignment(seq1, seq2);
		this.print_results();
		//System.out.println(this.out_seq1);System.out.println(this.out_seq2);

	}
	
	private void best_alignment(String seq1,String seq2) {
		long trace_start = System.nanoTime();
		//StringBuilder ultra efficient for working with strings 
		StringBuilder rev_seq1 = new StringBuilder(seq1); StringBuilder rev_seq2 = new StringBuilder(seq2);
		rev_seq1 = rev_seq1.reverse(); rev_seq2 = rev_seq2.reverse();
		int seq1_idx = 0; int seq2_idx = 0;
		
		//Start with the last position in the trace back matrix and go back from there.
		String how_align = trace_matrix[trace_matrix.length - 1][trace_matrix[0].length - 1];
		
		//use bytes because not as memory taxing
		byte[] out_seq1_arr = new byte[trace_matrix.length+trace_matrix[0].length]; byte[] out_seq2_arr = new byte[trace_matrix.length+trace_matrix[0].length];
		int align_idx = 0;
		int i = trace_matrix.length - 1;
		int j = trace_matrix[0].length - 1;
		
		while(how_align != "done") {
			if(how_align == "diag") {
				out_seq1_arr[align_idx] = (byte) rev_seq1.charAt(seq1_idx); out_seq2_arr[align_idx] = (byte) rev_seq2.charAt(seq2_idx);
				i--; seq1_idx++;
				j--; seq2_idx++;
			}
			else if(how_align == "up") {
				out_seq1_arr[align_idx] = (byte) rev_seq1.charAt(seq1_idx); out_seq2_arr[align_idx] = '-';
				i--; seq1_idx++;
			}
			else if(how_align == "left") {
				out_seq1_arr[align_idx] = '-'; out_seq2_arr[align_idx] = (byte) rev_seq2.charAt(seq2_idx);
				j--; seq2_idx++;
			}
			how_align = trace_matrix[i][j];
			align_idx++;
		}
		
		long best_align = System.nanoTime();
		System.out.println("Time to create alignment: "+(best_align - trace_start)/1000000000.0 + " seconds");
		
		String out_seq1 = trimmer(new String(out_seq1_arr)); String out_seq2 = trimmer(new String(out_seq2_arr)); 
		
		long trim_time = System.nanoTime();
		System.out.println("Time to trim strings: "+(trim_time - best_align)/1000000000.0 + " seconds");
		StringBuilder seq1_sb = new StringBuilder(out_seq1); StringBuilder seq2_sb = new StringBuilder(out_seq2);
		
		//Make my own trim method. should be quicker
		this.out_seq1 = seq1_sb.reverse().toString(); this.out_seq2 = seq2_sb.reverse().toString();
		long save_align = System.nanoTime();
		
		System.out.println("Time to save alignment: "+(save_align - best_align)/1000000000.0 + " seconds");
	}
	
	private static String trimmer(String str) {
		StringBuilder out = new StringBuilder("");
		int i = 0;
		while((int) str.charAt(i) != 0) {
			out.append(str.charAt(i));
			i++;
		}
		return out.toString();
	}
	
	private static int find_max(int[] arr) {
		int max = Integer.MIN_VALUE;
		for(int i:arr) {
			if(i > max) {
				max = i;
			}
		}
		
		return max;
	}
	
	public void print_score_mat() {
		for(int i = 0;i<score_matrix.length;i++) {
			for(int j = 0;j<score_matrix[0].length;j++) {
				System.out.print(score_matrix[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	public void print_trace_mat() {
		for(int i = 0;i<score_matrix.length;i++) {
			for(int j = 0;j<score_matrix[0].length;j++) {
				System.out.print(trace_matrix[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	public void print_results() {
		StringBuilder seq1 = new StringBuilder(this.out_seq1); StringBuilder seq2 = new StringBuilder(this.out_seq2);
		StringBuilder buffer = new StringBuilder("");
		int comp = 0;
		for(int i = 0; i< seq1.length();i++) {			
			comp = seq1.charAt(i) == seq2.charAt(i) ? 1:0;
				switch(comp) {
				case 1:
					buffer.append('*');
					break;
				case 0:
					buffer.append(' ');
					break;
				}
			}
		int j = 0;
		
		
		while(j<seq1.length()) {
			if(j%50 == 0) {
				seq1.insert(j,',');
				seq2.insert(j,',');
				buffer.insert(j,',');
			}
			j++;
		}
		
		
		
		String[] seq1_str = seq1.toString().split(",");
		String[] seq2_str = seq2.toString().split(",");
		String[] buff_str = buffer.toString().split(",");
		
	
		
		for(int i = 0; i < seq1_str.length;i++) {
			System.out.println(seq1_str[i]);
			System.out.println(buff_str[i]);
			System.out.println(seq2_str[i]);
			System.out.println();
		}
		
	}
}

class InputThread extends Thread{
	BufferedReader buff;
	String fin;
	
	public BufferedReader getBuff() {
		return buff;
	}

	public String getFin() {
		return fin;
	}

	public InputThread(BufferedReader buff) throws IOException {
		this.buff = buff;
	}
	
	public void run() {
		  String seq_line = "";
		try {
			seq_line = buff.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		  String seq_final = "";
		  
		  while (!(seq_line.isEmpty())) {
			  
			    if(seq_line.charAt(0) == '>') {
	            	try {
						seq_line = buff.readLine();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
	            	continue;
	            }
			  
	            else if(seq_line.charAt(0) != '>'){
	            	seq_final = seq_final + seq_line;
            		try {
						seq_line = buff.readLine();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
	            }
		  }
		  System.out.println("Thread "+this.getId()+" is done!");
		  fin = seq_final;
	}
	
}
