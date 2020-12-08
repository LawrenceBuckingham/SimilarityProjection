#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <algorithm>

#include "Array.hpp"

using namespace std;

int main(int argc, char **argv_) {

	// prevent accidental use of printf to emit results.
	int printf = 0;

	if ( argc < 5 ) {
		fprintf(stderr, "Usage: %s qrelsFile rankingFile summaryFile ignoreMissing=true|false [interpolated_precision_points=11]\n", argv_[0]);
		exit(1);
	}

	const char * qrelsFile = argv_[1];
	const char * rankingFile = argv_[2];
	const char * summaryFile = argv_[3];

	bool ignoreMissing = strcmp(argv_[4],"true") == 0;

	int iprecPoints = 11;

	if ( argc >= 6 ) {
		int numConverted = sscanf( argv_[5], "%d", &iprecPoints);

		if ( numConverted != 1 ) {
			fprintf(stderr, "Qrels file %s does not exist.\n", qrelsFile);
			exit(1);
		}
	}

	// fprintf(stderr, "Reading qrels\n");
	
	FILE *fp;
	fp = fopen(qrelsFile, "r");

	if ( !fp ) {
		fprintf(stderr, "Qrels file %s does not exist.\n", qrelsFile);
		exit(1);
	}

	size_t topicCount = 0, docCount = 0;
	size_t num_rel = 0;
	
	unordered_map<string, size_t> topicIds;
	vector<string> topicNames;

	unordered_map<string, size_t> docIds;
	vector<string> docNames;
	
	vector<unordered_set<size_t>> qrels;
	vector<size_t> per_topic_rel;

	for ( ;; ) {
		char topic_[256], doc_[256];
		int relevant;
		int items = fscanf(fp, "%s %*s %s %d\n", topic_, doc_, &relevant);
		if ( items < 3 ) {
			// fprintf(stderr, "%d items read. Finishing.\n", items);
			break;
		}

		string topic(topic_);
		string doc(doc_);
		size_t topicId;
		size_t docId;

		auto topicPos = topicIds.find(topic);

		if ( topicPos == topicIds.end() ) {
			topicId = topicCount++;
			topicIds.emplace(topic, topicId);
			topicNames.push_back(topic);
			per_topic_rel.push_back(0);
			qrels.resize(qrels.size() + 1);
		}
		else {
			topicId = topicPos->second;
		}

		auto docPos = docIds.find(doc);

		if ( docPos == docIds.end() ) {
			docId = docCount++;
			docIds.emplace(doc, docId);
			docNames.push_back(doc);
		}
		else {
			docId = docPos->second;
		}

		if ( relevant >= 1 ) {
			qrels[topicId].insert(docId);
			per_topic_rel[topicId]++;
			num_rel++;
		}
	}
	fclose(fp);
	// fprintf(stderr, "topicCount: %u\n", (unsigned) topicCount);

	// fprintf(stderr, "Reading rankings...\n");
	fp = fopen(rankingFile, "r");

	if ( !fp ) {
		fprintf(stderr, "Rankings file %s does not exist.\n", rankingFile);
		exit(1);
	}

	vector<bool> retrievedResultsFor(topicCount + 1);

	vector<vector<tuple<double, double, size_t>>> trec(topicCount + 1); // score, rank, docId
	vector<size_t> num_ret(topicCount + 1, 0);
	vector<size_t> num_rel_ret(topicCount + 1, 0);

	for ( ;; ) {
		char topic_[256], doc_[256];
		double rank, score;
		if ( fscanf(fp, "%s %*s %s %lf %lf %*[^\n]\n", topic_, doc_, &rank, &score) < 4 ) break;

		string topic(topic_);
		string doc(doc_);

		if ( topicIds.find(topic) != topicIds.end() ) {
			size_t topicId = topicIds[topic];

			size_t docId;
			auto docPos = docIds.find(doc);

			if ( docPos == docIds.end() ) {
				docId = docCount++;
				docIds.emplace(doc, docId);
				docNames.push_back(doc);
			}
			else {
				docId = docPos->second;
			}

			trec[topicId].push_back(make_tuple(-score, rank, docId));
			retrievedResultsFor[topicId] = true;

			// Update returned document count for individual query and total.
			num_ret[topicId] ++;
			num_ret[topicCount] ++;

			if ( qrels[topicId].count(docId) ) {
				// Update relevant document count for individual query and total.
				num_rel_ret[topicId] ++;
				num_rel_ret[topicCount] ++;
			}
		}
	}

	//for (auto t : trec[topicId]) {
	//  fprintf(stderr, "~ %s 0 %s %f %f\n", topicNames[topicId].c_str(), docNames[get<2>(t)].c_str(), get<1>(t), 0 - get<0>(t));
	//}


	fclose(fp);
	// fprintf(stderr, "Computing interpolated precision\n");
	vector<double> averagePrecisions(topicCount + 1);
	vector<vector<double>> iprec(topicCount);

#pragma omp parallel for
	for ( int i = 0; i < (int) topicCount; i++ ) {
		
		if ( (! retrievedResultsFor[i]) && ignoreMissing ) continue; 
		 
		sort(begin(trec[i]), end(trec[i]));
		vector<bool> relList;
		for ( size_t j = 0; j < trec[i].size(); j++ ) {
			const auto &r = trec[i][j];
			relList.push_back(qrels[i].count(get<2>(r)) != 0);
		}
		size_t relCount = per_topic_rel[i];
		size_t relFound = 0;
		size_t totalFound = 0;
		map<double, double> rp;
		double averagePrecision = 0.0;
		for ( bool r : relList ) {
			if ( r ) relFound++;
			double recall = 1.0 * relFound / relCount;
			totalFound++;
			double precision = 1.0 * relFound / totalFound;
			if ( r ) {

				averagePrecision += precision;
			}
			rp[recall] = max(rp[recall], precision);
		}
		if ( relCount > 0 ) {
			averagePrecision /= relCount;
		}
		map<double, double> irp;
		for ( const auto &r : rp ) {
			double currentP = 0.0;
			for ( const auto &r2 : rp ) {
				if ( r2.first >= r.first ) {
					currentP = max(currentP, r2.second);
				}
			}
			irp[r.first] = currentP;
		}

		averagePrecisions[i] = averagePrecision;

		for ( size_t j = 0; j < iprecPoints; j++ ) {
			double recall = 1.0 * j / (iprecPoints - 1);
			double currentP = 0.0;
			for ( const auto &r2 : rp ) {
				if ( r2.first >= recall ) {
					currentP = max(currentP, r2.second);
				}
			}
			iprec[i].push_back(currentP);
		}
	}

	size_t topicRetCount = 0;

	// fprintf(stderr, "Averaging\n");
	double meanAveragePrecision = 0.0;
	vector<double> averageIprec(iprecPoints, 0.0);
	for ( size_t i = 0; i < topicCount; i++ ) {

		if ( (!retrievedResultsFor[i]) && ignoreMissing ) continue;

		meanAveragePrecision += averagePrecisions[i];
		topicRetCount ++;

		for ( size_t j = 0; j < iprecPoints; j++ ) {
			averageIprec[j] += iprec[i][j];
		}
	}

	meanAveragePrecision /= topicRetCount;

	for ( size_t j = 0; j < iprecPoints; j++ ) {
		averageIprec[j] /= topicRetCount;
	}

	FILE *summaryStream = fopen(summaryFile, "w" "b");

	if ( !summaryStream ) {
		fprintf(stderr, "Unable to open summary file '%s'\n", summaryFile);
		exit(1);
	}

	// Emit results for each individual query.
	for ( size_t i = 0; i < topicCount; i++ ) {

		if ( (!retrievedResultsFor[i]) && ignoreMissing ) continue;

		fprintf(summaryStream, "%s\t%s\t%u\n", "num_ret", topicNames[i].c_str(), (unsigned) num_ret[i]);
		fprintf(summaryStream, "%s\t%s\t%u\n", "num_rel", topicNames[i].c_str(), (unsigned) per_topic_rel[i]);
		fprintf(summaryStream, "%s\t%s\t%u\n", "num_rel_ret", topicNames[i].c_str(), (unsigned) num_rel_ret[i]);
		fprintf(summaryStream, "%s\t%s\t%.4f\n", "map", topicNames[i].c_str(), averagePrecisions[i]);

		for ( size_t j = 0; j < iprecPoints; j++ ) {
			double recall = 1.0 * j / (iprecPoints - 1);
			fprintf(summaryStream, "iprec_at_recall_%.2f\t%s\t%.4f\n", recall, topicNames[i].c_str(), iprec[i][j]);
		}
	}

	// Emit results for total.

	fprintf(summaryStream, "%s\tall\t%d\n", "num_q", (unsigned) topicCount);
	fprintf(summaryStream, "%s\tall\t%u\n", "num_ret", (unsigned) num_ret[topicCount]);
	fprintf(summaryStream, "%s\tall\t%u\n", "num_rel", (unsigned) num_rel);
	fprintf(summaryStream, "%s\tall\t%u\n", "num_rel_ret", (unsigned) num_rel_ret[topicCount]);
	fprintf(summaryStream, "%s\tall\t%.4f\n", "map", meanAveragePrecision);

	for ( size_t j = 0; j < iprecPoints; j++ ) {
		double recall = 1.0 * j / (iprecPoints - 1);
		fprintf(summaryStream, "iprec_at_recall_%.2f\tall\t%.4f\n", recall, averageIprec[j]);
	}

	fclose(summaryStream);

	// fprintf(stderr, "Finished.\n");
	return 0;
}

