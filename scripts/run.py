#!/usr/bin/env python

import sys, os, re, shutil, subprocess, snp_phylo_utils, random
from contextlib import contextmanager

@contextmanager
def cd(newdir):
	prevdir = os.getcwd()
	os.chdir(os.path.expanduser(newdir))
	try:
		yield
	finally:
		os.chdir(prevdir)

def remove_results_without_reads(base_path, reads_dir):
	analr = os.path.join(base_path, "analysis_results")
	for x in os.listdir(analr):
		abspath = os.path.join(analr, x)
		if os.path.isdir(abspath) == True:
			matches = len([y for y in os.listdir(reads_dir) if re.match(x, y)])
			if matches == 1:
				continue
			elif matches > 1:
				print("Something strange going on here. More than one read file for results directory {}. You should check this.".format(x))
			else:
				print("{} has no associated read file and will be deleted.".format(x))
				shutil.rmtree(abspath)

def check_ref(old_reference):
	if len(os.listdir(old_reference)) > 1:
		print("Please make sure there is only one reference file in raw_data/reference_sequence directory.")
		sys.exit(1)
	new_reference = os.path.join(old_reference, os.listdir(old_reference)[0])
	if os.path.splitext(new_reference)[1] == ".gbk" or os.path.splitext(new_reference)[1] == ".gb" or os.path.splitext(new_reference)[1] == ".gbff":
		return new_reference, ".gbk"
	elif os.path.splitext(new_reference)[1] == ".fasta" or os.path.splitext(new_reference)[1] == ".fa" or os.path.splitext(new_reference)[1] == ".fna":
		return new_reference, ".fasta"
	else:
		print("Can't id reference sequence format. Ensure you have a Genbank (.gb, .gbk, .gbff) or Fasta (.fasta,.fa,.fna) formatted sequence.")
		sys.exit(1)

def remove_degenerate_bases(reference, refsuff, outdir, threads):
	if os.path.isfile(os.path.join(outdir, "reference{}".format(refsuff))) == False:
		snp_phylo_utils.remove_degenerate_bases(reference, os.path.join(outdir, "reference"), threads)
		# subprocess.run(["python", "remove_degenerate_bases.py", reference, os.path.join(outdir, "reference"), threads], check = True, stderr=subprocess.STDOUT)
	else:
		print("Reference sequence has already been processed. Skipping this step...")


def run_snippy(read_path, threads, outdir, reference):
	name = snp_phylo_utils.get_samplename(read_path)
	read_format = snp_phylo_utils.id_read_format(read_path)
	if read_format == 'fastq':
		try:
			subprocess.run(['snippy', '--cpus', str(threads), '--prefix', name, '--outdir', os.path.join(outdir, name), '--ref', reference, '--peil', read_path], check = True, stderr=subprocess.STDOUT)
		except Exception as e:
			print(e)
			print("Error running snippy on fastq reads file: {}. Please check your files and error output.".format(os.path.basename(read_path)))
			raise
	else:
		try:
			subprocess.run(['snippy', '--cpus', str(threads), '--prefix', name, '--outdir', os.path.join(outdir, name), '--ref', reference, '--ctgs', read_path], check = True, stderr=subprocess.STDOUT)
		except Exception as e:
			print(e)
			print("Error running snippy on fasta assembled genome file: {}. Please check your files and error output.".format(os.path.basename(read_path)))
			raise


def hpc_run(read_path, outdir, reference):
	run_snippy(read_path, 1, outdir, reference)
	sys.exit(0)


def run_snippy_core(outdir, id_threshold, sample_names):
	if os.path.isfile(os.path.join(outdir, "core.full.aln")) == False:
		with cd(outdir):
			subprocess.run((["snippy-core"]+sample_names), check = True, stderr=subprocess.STDOUT)
			os.makedirs("../excluded_sequences/poor_ref_alignment", exist_ok=True)
			bad_seqs = snp_phylo_utils.ditch_distant_sequences('core.txt', id_threshold)
			subprocess.run(["snippy-core"]+[x for x in sample_names if (x in bad_seqs) == False], check = True, stderr=subprocess.STDOUT)
	else:
		print("Core alignment has already been generated. Skipping this step...")

def Ns2gaps(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "core.full.trimmed.aln")) == False:
		with cd(outdir):
			snp_phylo_utils.replace_Ns_with_gaps("core.full.aln", threads, "core.full.trimmed.aln")
	else:
		print("Core alignment has already been trimmed. Skipping this step...")

def count_total_snps(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "total.snp_counts.txt")) == False:
		with cd(outdir):
			snp_phylo_utils.count_snps('core.full.trimmed.aln', 'total.snp_counts.txt', 'total.snp_counts.stats', threads)
	else:
		print("Total snp count matrix has already been generated. Skipping this step...")

def filter_recombinant_positions(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "filtered_core_aln.final_tree.tre")) == False:
		with cd(outdir):
			subprocess.run(["bash","../scripts/filter_recomb.sh",str(threads)], check = True, stderr=subprocess.STDOUT)
	else:
		print("Recombination filtering by gubbins has already been done. Skipping this step...")

def visualise_recombination_predictions(outdir):
	if os.path.isfile(os.path.join(outdir, "filtered_core_aln.recombination_predictions.pdf")) == False:
		with cd(outdir):
			subprocess.run(["bash","../scripts/visualise_recomb.sh"], check = True, stderr=subprocess.STDOUT)
	else:
		print("Recombination prediction visualisation has already been generated. Skipping this step...")

def count_core_snps(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "filtered_core_aln.gubbins_filtered.snp_counts.txt")) == False:
		with cd(outdir):
			snp_phylo_utils.count_snps('filtered_core_aln.filtered_polymorphic_sites.fasta', 'filtered_core_aln.gubbins_filtered.snp_counts.txt', 'filtered_core_aln.gubbins_filtered.snp_counts.stats', threads)
	else:
		print("Core snp count matrix has already been generated. Skipping this step...")

def bootstrap_tree(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "RAxML_bootstrap.filtered_core_aln.bootstrap")) == False:
		with cd(outdir):
			subprocess.run(["raxmlHPC-PTHREADS","-T",str(threads),"-m","GTRCAT","-V","-p",str(random.randint(10000,99999)),"-b",str(random.randint(10000,99999)),"-#","100","-s","filtered_core_aln.filtered_polymorphic_sites.fasta","-n","filtered_core_aln.bootstrap"], check = True, stderr=subprocess.STDOUT)
	else:
		print("Bootstrap RAxML trees have already been generated. Skipping this step...")

def gen_final_tree(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "RAxML_bipartitionsBranchLabels.filtered_core_aln.final")) == False:
		with cd(outdir):
			subprocess.run(["raxmlHPC-PTHREADS","-T",str(threads),"-m","GTRCAT","-V","-p",str(random.randint(10000,99999)),"-f","b","-t","filtered_core_aln.final_tree.tre","-z","RAxML_bootstrap.filtered_core_aln.bootstrap","-n","filtered_core_aln.final"], check = True, stderr=subprocess.STDOUT)
	else:
		print("Final RAxML tree has already been generated. Skipping this step...")

def gen_final_rooted_tree(outdir, threads):
	if os.path.isfile(os.path.join(outdir, "RAxML_rootedTree.filtered_core_aln.finalrooted")) == False:
		with cd(outdir):
			subprocess.run(["raxmlHPC-PTHREADS","-T",str(threads),"-m","GTRCAT","-V","-f","I","-t","RAxML_bipartitionsBranchLabels.filtered_core_aln.final","-n","filtered_core_aln.finalrooted"], check = True, stderr=subprocess.STDOUT)
	else:
		print("Final RAxML rooted tree has already been generated. Skipping this step...")

def cleanup(outdir):
	with cd(outdir):
		os.mkdir("intermediate_files")
		i_files = ["core.full.trimmed.aln",
				"core.full.trimmed.aln.seq.joint.txt",
				"core.branch_base_reconstruction.embl",
				"core.filtered_polymorphic_sites.fasta",
				"core.filtered_polymorphic_sites.phylip",
				"core.final_tree.tre",
				"core.per_branch_statistics.csv",
				"core.gubbins_filtered.trimmed.aln",
				"core.gubbins_filtered.trimmed.aln.reduced",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.0",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.1",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.2",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.3",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.4",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.5",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.0",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.1",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.2",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.3",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.4",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.5",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.6",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.6",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.7",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.8",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.9",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.10",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.11",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.12",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.7",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.8",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.9",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.10",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.11",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.12",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.13",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.13",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.14",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.15",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.16",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.17",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.18",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.14",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.15",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.16",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.17",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.18",
				"RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.19",
				"RAxML_bestTree.core.gubbins_filtered.trimmed",
				"RAxML_info.core.gubbins_filtered.trimmed",
				"RAxML_log.core.gubbins_filtered.trimmed.RUN.19",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.0",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.1",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.2",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.3",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.4",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.5",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.6",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.7",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.8",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.9",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.10",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.11",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.12",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.13",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.14",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.15",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.16",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.17",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.18",
				"RAxML_result.core.gubbins_filtered.trimmed.RUN.19",
				"RAxML_bootstrap.core.gubbins_filtered.trimmed.bootstrap",
				"RAxML_info.core.gubbins_filtered.trimmed.bootstrap",
				"RAxML_bipartitions.core.gubbins_filtered.trimmed.final",
				"RAxML_bipartitionsBranchLabels.core.gubbins_filtered.trimmed.final",
				"RAxML_info.core.gubbins_filtered.trimmed.final",
				"RAxML_info.core.gubbins_filtered.trimmed.finalrooted"
				]
		for file in i_files:
			if os.path.exists(file):
				shutil.move(file, "intermediate_files")
		


def main():
	base_path = os.path.abspath(sys.argv[1])
	reads_dir = os.path.abspath(sys.argv[2])
	old_reference = os.path.abspath(sys.argv[3])
	outdir = os.path.abspath(sys.argv[4])
	threads = int(sys.argv[5])
	hpc = sys.argv[6]
	id_threshold = sys.argv[7]
	# for x in sys.argv:
		# print(x)
	if sys.argv[8] != 'all':
		hpc_run(sys.argv[8], outdir, old_reference)
	temp_dir = os.path.join(base_path, '.temp')
	remove_results_without_reads(base_path, reads_dir)
	reference, refsuff = check_ref(old_reference)
	remove_degenerate_bases(reference, refsuff, outdir, threads)
	reference = os.path.join(outdir, "reference{}".format(refsuff))
	sample_names = [ snp_phylo_utils.get_samplename(path) for path in os.listdir(reads_dir) if os.path.isfile(os.path.join(reads_dir, path)) ]
	print(sample_names)
	snippy_not_done = [ os.path.join(reads_dir, path) for path in os.listdir(reads_dir) if (snp_phylo_utils.get_samplename(path) in os.listdir(outdir)) == False ]
	if len(snippy_not_done) != 0:
		if hpc == 'true':
			mem = sys.argv[9]
			alignment_ids = []
			for i, read_path in enumerate(snippy_not_done):
				job_id = snp_phylo_utils.submit_hpc_script(temp_dir, read_path, outdir, reference, '1', '4gb', 'smallq', base_path, id_threshold)
				alignment_ids.append(job_id)
			print(alignment_ids)
			snp_phylo_utils.submit_tree_hpc_script(temp_dir, read_path, outdir, old_reference, threads, mem, 'workq', base_path, id_threshold, alignment_ids)
			# print("Alignment stage is now running as {} individual submitted jobs. Check the status of your jobs with the qstat command. \n\nWhen all of your jobs are completed (no jobs shown on qstat) they please restart this pipeline with the following command:\n\n{}\n\nto complete the analysis.".format((i+1),sys.argv[0]))
			sys.exit(0)
		else:
			for read_path in snippy_not_done:
				run_snippy(read_path, threads, outdir, reference)
	run_snippy_core(outdir, id_threshold, sample_names)
	Ns2gaps(outdir, threads)
	count_total_snps(outdir, threads)
	filter_recombinant_positions(outdir, threads)
	visualise_recombination_predictions(outdir)
	count_core_snps(outdir, threads)
	bootstrap_tree(outdir, threads)
	gen_final_tree(outdir, threads)
	gen_final_rooted_tree(outdir, threads)
	cleanup(outdir)

if __name__ == '__main__':
	main()
