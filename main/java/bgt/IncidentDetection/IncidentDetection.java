package bgt.IncidentDetection;

import bgt.IncidentDetection.Models.EvaluateResult;
import bgt.IncidentDetection.Models.LabeledRecord;
import bgt.IncidentDetection.Models.SegDistribution;
import org.apache.commons.io.FileUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by admin on 2016/12/2.
 */
public class IncidentDetection {
    public List<SegDistribution> segDistributionList;
    public HashMap<String, SegDistribution> segDistributionMap;
    public static int TIME_WINDOW = 6; //can be 60, 300, 600, 1800, 3600

    public void createDir(String dirName){
    	File theDir = new File(dirName);
    	// if the directory does not exist, create it
    	if (!theDir.exists()) {
    	    System.out.println("creating directory: " + theDir.getName());
    	    boolean result = false;

    	    try{
    	        theDir.mkdir();
    	        result = true;
    	    } 
    	    catch(SecurityException se){
    	        //handle it
    	    }        
    	    if(result) {    
    	        System.out.println("DIR created");  
    	    }
    	}
    }
    public void parseDistribution(String distribution_filename) throws FileNotFoundException, IOException{
        segDistributionList = new ArrayList<>();
        segDistributionMap = new HashMap<>();
        InputStreamReader is = new InputStreamReader(new FileInputStream(distribution_filename));
        BufferedReader br = new BufferedReader(is);
        String line = null;
        int line_count = 1;
        while ((line = br.readLine()) != null){
            if (line_count == 1){   // info about # of mixtures
                String[] line_split = line.split(" ");
                SegDistribution.NUM_MIXTURE = Integer.valueOf(line_split[line_split.length-1]);
                line_count++;
            }else if (line_count == 3){
                SegDistribution.lambdas = new double[SegDistribution.NUM_MIXTURE];
                SegDistribution.deltas = new double[SegDistribution.NUM_MIXTURE][SegDistribution.NUM_MIXTURE];
                for (int i=0; i<SegDistribution.NUM_MIXTURE; i++){
                    line = br.readLine();
                    SegDistribution.lambdas[i] = Double.valueOf(line.split(",")[0]);
                    //System.out.println(SegDistribution.lambdas[i]);
                    line_count++;
                }
                for (int i=0; i<SegDistribution.NUM_MIXTURE; i++){
                    for (int j=0; j<SegDistribution.NUM_MIXTURE; j++){
                        SegDistribution.deltas[i][j] = KL_Distance(SegDistribution.lambdas[i], SegDistribution.lambdas[j]);
                        //System.out.println(SegDistribution.deltas[i][j]);
                    }
                }// DELTAS ON SIIS MINGI KL_DISTANCE BETWEEN THE K DISTRIBUTION.
            }else if(line_count >= SegDistribution.NUM_MIXTURE+6){
                String[] line_split = line.split(",");
                String seg_id = line_split[0];
                double[] pis = new double[SegDistribution.NUM_MIXTURE];
                for (int i=2; i<2+SegDistribution.NUM_MIXTURE; i++){
                    pis[i-2] = Double.valueOf(line_split[i]);
                    //System.out.println(i+","+pis[i-2]);
                }
                double[] sigmas = new double[SegDistribution.NUM_MIXTURE];
                for (int i=0; i< SegDistribution.NUM_MIXTURE; i++){
                    double sigma = 0.0;
                    for (int j=0; j < SegDistribution.NUM_MIXTURE; j++){
                        sigma += pis[j]*SegDistribution.deltas[j][i];
                    }
                    sigmas[i] = sigma;
                }
                SegDistribution segDistribution = new SegDistribution(seg_id, pis, sigmas);
                segDistributionList.add(segDistribution);
                segDistributionMap.put(seg_id, segDistribution);
                line_count++;
            }else{
                line_count++;
            }
        }
        br.close();
        is.close();
    }

    public double factorial(int x) {
        double fact = 1;
        for (int i = 2; i <= x; i++) {
            fact *= i;
        }
        return fact;
    }
    
    public double poissonDistribution(double lambda, int x){
        return (Math.pow(lambda, x)*Math.exp(-lambda))/factorial(x);
    }

    public double KL_Distance(double lambda_k, double lambda_l){
        return lambda_l - lambda_k + lambda_k*Math.log(lambda_k/lambda_l);
    }

    public double evalDivergence(String seg_id, int speed){
        if (segDistributionMap.containsKey(seg_id)){
            SegDistribution segDistribution = segDistributionMap.get(seg_id);
            // look for the Usual Traffic State
            int delta_s = 0;
            double max_pi_s = segDistribution.pis[0];
            //System.out.println(max_pi_s);
            for (int i=1; i<SegDistribution.NUM_MIXTURE; i++){
                if (segDistribution.pis[i] > max_pi_s){
                    delta_s = i;
                    max_pi_s = segDistribution.pis[i];
                }
            }

            // look for the Current Traffic State
            int delta_sx = 0;
            double max_pi_sx = segDistribution.pis[0] * poissonDistribution(SegDistribution.lambdas[0], speed);
            for (int i=1; i<SegDistribution.NUM_MIXTURE; i++){
                double this_pi_sx = segDistribution.pis[i]* poissonDistribution(SegDistribution.lambdas[i], speed);
                if (this_pi_sx > max_pi_sx){
                    delta_sx = i;
                    max_pi_sx = this_pi_sx;
                }
            }
            double lambda_s = SegDistribution.lambdas[delta_s];
            double lambda_sx = SegDistribution.lambdas[delta_sx];
            //System.out.println(delta_s+","+delta_sx);
            return KL_Distance(lambda_s, lambda_sx);
        }
        return Double.MAX_VALUE;
    }

    public double evalDivergence_weighted(String seg_id, int speed){
        if (segDistributionMap.containsKey(seg_id)){
            SegDistribution segDistribution = segDistributionMap.get(seg_id);
            double[] sigmas = segDistribution.sigmas;
            double[] pis = new double[SegDistribution.NUM_MIXTURE];
            double pi_sum = 0.0;
            for (int i=0; i<SegDistribution.NUM_MIXTURE; i++){
                pis[i] = segDistribution.pis[i]*poissonDistribution(SegDistribution.lambdas[i], speed);
                pi_sum += pis[i];
            }
            double result = 0.0;
            for (int i=0; i<SegDistribution.NUM_MIXTURE; i++){
                result += sigmas[i]*(pis[i]/pi_sum);
            }
            return result;
        }
        return Double.MAX_VALUE;
    }

    /*public List<EvaluateResult> detectOnStream(String seg_id, List<LabeledRecord> recList){
        // return list of anomaly values
        List<EvaluateResult> streamResult = new ArrayList<>();
        int totalRecNum = recList.size();
        if (totalRecNum == 1){
            double divergence = evalDivergence(seg_id, recList.get(0).speed);
            double divergence_weighted = evalDivergence_weighted(seg_id, recList.get(0).speed);
            EvaluateResult evalRes = new EvaluateResult(seg_id, divergence, divergence_weighted, recList.get(0).acc_flag);
            streamResult.add(evalRes);
        }else{
            int windowStart = 0;
            while (windowStart < recList.size()){   // until start from the last record
                int windowEnd = windowStart;
                
                while ( windowEnd+1 < recList.size() &&
                        ((recList.get(windowEnd+1).time - recList.get(windowStart).time)*60) < TIME_WINDOW){
                    windowEnd++;
                    //System.out.println((recList.get(windowEnd+1).time +","+ recList.get(windowStart).time+","+(recList.get(windowEnd+1).time - recList.get(windowStart).time)*60.0));
                    //System.out.println("window end: "+ recList.get(windowEnd+1).time + " window start: "+ recList.get(windowStart).time);
                }
                //System.out.println("window end: "+ windowEnd+ " now: "+ windowStart);

                double divSum = 0;
                double divSum_weighted = 0;
                int acc_flag = 0;
                for (int i=windowStart; i<=windowEnd; i++){
                    divSum += evalDivergence(seg_id, recList.get(i).speed);
                    divSum_weighted += evalDivergence_weighted(seg_id, recList.get(i).speed);
                    if (recList.get(i).acc_flag >= 1) acc_flag = 1;
                    //if(divSum==0&&recList.get(i).acc_flag>0){
                    	//System.out.println("Different State "+seg_id+" Time: "+ recList.get(i).time+", Speed: "+recList.get(i).speed+", Acc: "+ recList.get(i).acc_flag);
                    //}
                }
                int windowLen = windowEnd-windowStart+1;
                //System.out.println("window len: "+ windowLen);
                EvaluateResult evalRes = new EvaluateResult(seg_id, divSum/windowLen, divSum_weighted/windowLen, acc_flag);
                
                streamResult.add(evalRes);

                windowStart++;  // next window
                //System.out.println(""+evalRes.anomaly_val+","+ evalRes.anomaly_val_weighted);
            }
            //System.out.println("detectOnstream: all loop finished");
        }
        return streamResult;
    }*/
    

    /*public void parseLabeledRecords(String label_filename_prefix, int fold) throws FileNotFoundException, IOException{
        List<Long> acc_route_list = new ArrayList<>();   // List of route numbers that has at least one positive record
        FileUtils.cleanDirectory(new File("data/incident_detection/TW_"+TIME_WINDOW+"/fold_"+fold+"/"));

        InputStreamReader is = new InputStreamReader(new FileInputStream(label_filename_prefix+"_list.txt"));
        BufferedReader br = new BufferedReader(is);
        FileWriter fw;
        String line = null;
        while ((line = br.readLine()) != null){
            String[] line_split = line.split("\t");
            acc_route_list.add(Long.valueOf(line_split[0])); //
        }
        br.close();
        is.close();
        for (long route_num : acc_route_list){// Reads in way files
            is = new InputStreamReader(new FileInputStream(label_filename_prefix+"_"+route_num+".txt"));
            br = new BufferedReader(is);
            List<EvaluateResult> routeResult = new ArrayList<>();
            String seg_id = null;
            List<LabeledRecord> recList = new ArrayList<>();   // one day's record stream for one segment
            while ((line = br.readLine()) != null){
                if (!line.equals("")){  // not empty line
                    String[] line_split = line.split("\t");
                    if (seg_id == null) seg_id = line_split[0];
                    LabeledRecord lRec = new LabeledRecord(Double.valueOf(line_split[1]), Integer.valueOf(line_split[2]), Integer.valueOf(line_split[3]));
                    recList.add(lRec);
                }else{
                    routeResult.addAll(detectOnStream(seg_id, recList));
                    seg_id = null;
                    recList.clear();
                }
            }
            routeResult.addAll(detectOnStream(seg_id, recList));    // last stream
            br.close();
            is.close();

            fw = new FileWriter("data/incident_detection/TW_"+TIME_WINDOW+"/fold_"+fold+"/detect_result_"+route_num+".txt",true);
            double min_anomaly_value = Double.MAX_VALUE, min_anomaly_value_weighted = Double.MAX_VALUE;
            double max_anomaly_value = Double.MIN_VALUE, max_anomaly_value_weighted = Double.MIN_VALUE;
            for (EvaluateResult evalRes : routeResult){
                min_anomaly_value = evalRes.anomaly_val < min_anomaly_value ? evalRes.anomaly_val : min_anomaly_value;
                min_anomaly_value_weighted = evalRes.anomaly_val_weighted < min_anomaly_value_weighted ?
                        evalRes.anomaly_val_weighted : min_anomaly_value_weighted;
                max_anomaly_value = (evalRes.anomaly_val > max_anomaly_value && evalRes.anomaly_val!=Double.MAX_VALUE && Double.isFinite(evalRes.anomaly_val))? evalRes.anomaly_val : max_anomaly_value;
                max_anomaly_value_weighted = (evalRes.anomaly_val_weighted > max_anomaly_value_weighted && evalRes.anomaly_val_weighted!=Double.MAX_VALUE && Double.isFinite(evalRes.anomaly_val_weighted))?
                        evalRes.anomaly_val_weighted : max_anomaly_value_weighted;
            }
            int acc_num = 0,total_num = 0;
            for (EvaluateResult evalRes : routeResult){
                if (evalRes.acc_label == 1){
                    acc_num++;
                }
                total_num++;
            }
            fw.write("DR\tFAR\n");
            fw.write("KL Result\n");
            double threshold = min_anomaly_value;
            double step = (max_anomaly_value - min_anomaly_value) / 1000.0 + Double.MIN_VALUE;
            int detected_num, false_alarm_num;
            while (threshold <= max_anomaly_value){
                detected_num = 0; false_alarm_num = 0;
                for (EvaluateResult evalRes : routeResult){
                    if (evalRes.acc_label == 1 && evalRes.anomaly_val >= threshold){
                        detected_num++;
                    }
                    if (evalRes.acc_label == 0 && evalRes.anomaly_val >= threshold){
                        false_alarm_num++;
                    }
                }
                
                double DR = (double) detected_num / (double) acc_num;
                double FAR = (double) false_alarm_num / (double) total_num;
                fw.write(DR+"\t"+FAR+"\n");
                threshold += step;
            }
            fw.write("Weighted KL Result\n");

            threshold = min_anomaly_value_weighted;
            step = (max_anomaly_value_weighted - min_anomaly_value_weighted) / 1000.0+ Double.MIN_VALUE;
            while (threshold <= max_anomaly_value_weighted){
                detected_num = 0; false_alarm_num = 0;
                for (EvaluateResult evalRes : routeResult){
                    if (evalRes.acc_label == 1 && evalRes.anomaly_val_weighted >= threshold){
                            detected_num++;
                    }
                    if (evalRes.acc_label == 0 && evalRes.anomaly_val_weighted >= threshold){
                        false_alarm_num++;
                    }
                }
                double DR = (double) detected_num / (double) acc_num;
                double FAR = (double) false_alarm_num / (double) total_num;
                fw.write(DR+"\t"+FAR+"\n");
                threshold += step;
            }
            fw.close();
            br.close();
            is.close();
        }
    }*/
    
    public void parseLabeledRecords(String label_filename_prefix, int k) throws FileNotFoundException, IOException{
        List<Long> acc_route_list = new ArrayList<>();   // List of route numbers that has at least one positive record
        this.createDir("data/incident_detection/TW_"+TIME_WINDOW+"/");
        if(k==2)
        	FileUtils.cleanDirectory(new File("data/incident_detection/TW_"+TIME_WINDOW+"/"));

        InputStreamReader is;
        BufferedReader br;
        FileWriter fw;
        String line = null;
        //for (long route_num : acc_route_list){// Reads in way files
            is = new InputStreamReader(new FileInputStream(label_filename_prefix));
            br = new BufferedReader(is);
            List<EvaluateResult> routeResult = new ArrayList<>();
            String seg_id = null;
            List<LabeledRecord> recList = new ArrayList<>();   // one day's record stream for one segment
            String prevDay = null;
            //String prevId = null;
            /*while ((line = br.readLine()) != null){
                if (!line.equals("")){  // not empty line
                    String[] line_split = line.split("\t");
                    if (seg_id == null) seg_id = line_split[0];
                    LabeledRecord lRec = new LabeledRecord(Double.valueOf(line_split[1]), Integer.valueOf(line_split[2]), Integer.valueOf(line_split[3]));
                    recList.add(lRec);
                }else{
                    routeResult.addAll(detectOnStream(seg_id, recList));
                    seg_id = null;
                    recList.clear();
                }
            }*/
            /*while ((line = br.readLine()) != null){
            	String[] line_split = line.split("\t");
            	if(seg_id == null){
            		//prevDay = line_split[4];
            		seg_id = line_split[0];
            	}
            	if (seg_id.equals(line_split[0])){  // not empty line
                    //seg_id = line_split[0];
                    //System.out.println(Double.valueOf(line_split[1])+","+Integer.valueOf(line_split[2])+","+Integer.valueOf(line_split[3]));
                    //LabeledRecord lRec = new LabeledRecord(Double.valueOf(line_split[1]), Integer.valueOf(line_split[2]), Integer.valueOf(line_split[3]));
                    //EvaluateResult evalRes = new EvaluateResult(seg_id, divSum/windowLen, divSum_weighted/windowLen, acc_flag);
                    //routeResult.add(arg0)
                    //if(Integer.valueOf(line_split[3])>0)
                    		//System.out.println("Label_record: "+Double.valueOf(line_split[1])+", "+Integer.valueOf(line_split[2])+", "+Integer.valueOf(line_split[3]));
                    //recList.add(lRec);
                }else{
                	//System.out.println(recList.size());
                	//System.out.println(line_split[4]+','+prevDay);
                    //routeResult.addAll(detectOnStream(seg_id, recList));
                    recList.clear();
                    seg_id = line_split[0];
                    LabeledRecord lRec = new LabeledRecord(Double.valueOf(line_split[1]), Integer.valueOf(line_split[2]), Integer.valueOf(line_split[3]));
                    //System.out.println("New Day");
                    //if(Integer.valueOf(line_split[3])>0)
                    	//System.out.println("Label_record: "+Double.valueOf(line_split[1])+", "+Integer.valueOf(line_split[2])+", "+Integer.valueOf(line_split[3]));
                    recList.add(lRec);
                }
            }*/
            
            int prevAcc = 0;
            String prevSegId = "WRONG";
            while((line = br.readLine()) != null){
            	String[] line_split = line.split("\t");
            	LabeledRecord lRec = new LabeledRecord(line_split[0], Integer.valueOf(line_split[4]),Double.valueOf(line_split[1]), Integer.valueOf(line_split[2]), Integer.valueOf(line_split[3]));
         
            	if(prevAcc == 0 && lRec.acc_flag == 1){
            		lRec.isStart = 1;
            	}
            	else if(!prevSegId.equals(lRec.seg_id.split("_")[0])  && lRec.acc_flag == 1){
            		prevSegId = lRec.seg_id.split("_")[0];
            		//System.out.println(prevSegId);
            		lRec.isStart = 1;
            	}
            	//else if(prevAcc == 1 && lRec.acc)
            	else
            		lRec.isStart = 0;
            	
            	//if(lRec.isStart == 1)
            		//System.out.println(lRec.seg_id+','+ lRec.day+','+lRec.time+','+ lRec.speed+','+ lRec.acc_flag +"," + lRec.isStart);
            	
            	prevAcc = lRec.acc_flag;
            	recList.add(lRec);
            }
            System.out.println("Record Read");
            int windowStart = 0;
            for(LabeledRecord record : recList){
                int totalRecNum = recList.size();
                //while (windowStart < recList.size()){   // until start from the last record
                int windowEnd = windowStart;
                while ( windowEnd+1 < recList.size() &&
                        recList.get(windowEnd).seg_id.split("_")[0].equals(recList.get(windowStart).seg_id.split("_")[0]) &&
                        (windowEnd - windowStart + 1) < TIME_WINDOW){
                    windowEnd++;
                    //System.out.println((recList.get(windowEnd+1).time +","+ recList.get(windowStart).time+","+(recList.get(windowEnd+1).time - recList.get(windowStart).time)*60.0));
                    //System.out.println("window end: "+ recList.get(windowEnd+1).time + " window start: "+ recList.get(windowStart).time);
               
                }
                int windowLen = windowEnd-windowStart+1;
                //System.out.println("window len: "+ windowLen);
                double divSum = 0;
                double divSum_weighted = 0;
                int acc_flag = 0;
                int isStart = 0;
                for (int i=windowStart; i<=windowEnd; i++){
                    divSum += evalDivergence(record.seg_id, recList.get(i).speed);
                    divSum_weighted += evalDivergence_weighted(record.seg_id, recList.get(i).speed);
                    
                    //if(divSum==0&&recList.get(i).acc_flag>0){
                    	//System.out.println("Different State "+seg_id+" Time: "+ recList.get(i).time+", Speed: "+recList.get(i).speed+", Acc: "+ recList.get(i).acc_flag);
                    //}
                }
                if (record.acc_flag == 1) acc_flag = 1;
                if (record.isStart == 1) isStart =1;
                
               
                EvaluateResult evalRes = new EvaluateResult(record.seg_id, divSum/windowLen, divSum_weighted/windowLen, acc_flag, isStart);
                
                routeResult.add(evalRes);
                windowStart++;  // next window
                    //System.out.println(""+evalRes.anomaly_val+","+ evalRes.anomaly_val_weighted);
                //}
                //System.out.println("detectOnstream: all loop finished");
            }
            //routeResult.addAll(detectOnStream(seg_id, recList));    // last stream
            System.out.println("1 loop finished, "+ routeResult.size());
            br.close();
            is.close();
            
            //fw = new FileWriter("data/incident_detection/TW_"+TIME_WINDOW+"/detect_result_k="+k+".txt",true);
            double min_anomaly_value = Double.MAX_VALUE, min_anomaly_value_weighted = Double.MAX_VALUE;
            double max_anomaly_value = Double.MIN_VALUE, max_anomaly_value_weighted = Double.MIN_VALUE;
            for (EvaluateResult evalRes : routeResult){
                min_anomaly_value = evalRes.anomaly_val < min_anomaly_value ? evalRes.anomaly_val : min_anomaly_value;
                min_anomaly_value_weighted = evalRes.anomaly_val_weighted < min_anomaly_value_weighted ?
                        evalRes.anomaly_val_weighted : min_anomaly_value_weighted;
                max_anomaly_value = (evalRes.anomaly_val > max_anomaly_value && evalRes.anomaly_val!=Double.MAX_VALUE && Double.isFinite(evalRes.anomaly_val))? evalRes.anomaly_val : max_anomaly_value;
                max_anomaly_value_weighted = (evalRes.anomaly_val_weighted > max_anomaly_value_weighted && evalRes.anomaly_val_weighted!=Double.MAX_VALUE && Double.isFinite(evalRes.anomaly_val_weighted))?
                        evalRes.anomaly_val_weighted : max_anomaly_value_weighted;
                //System.out.println(""+min_anomaly_value+','+min_anomaly_value_weighted+','+max_anomaly_value+','+max_anomaly_value_weighted);
            }
            System.out.println("2 loop finished");
            int[] totalTime = new int[24];
            int[] accTime = new int[24];
            int[] normalTime = new int[24];
            int acc_num = 0,total_num = 0, normal_num = 0;
            for (EvaluateResult evalRes : routeResult){
            	int evalResTime = Integer.parseInt(evalRes.seg_id.split("_")[1]);
                if (evalRes.isStart == 1){
                    acc_num++;
                    accTime[evalResTime]++;
                }
                if (evalRes.acc_label == 0){
                	normal_num ++;
                	normalTime[evalResTime]++;
                }
                total_num++;
                totalTime[evalResTime]++;
            }
            System.out.println("incident num:"+acc_num);
            System.out.println("3 loop finished");
            /*fw.write("FARtotal\tDRtotal\t");
            for(int i=0;i<24;i++){
            	fw.write("FAR"+i+"\t"+"DR"+i+"\t");
            }
            fw.write("\n");
            //fw.write("KL Result\n");
            double threshold = min_anomaly_value;
            //double step = Double.MIN_VALUE;
            double step = (max_anomaly_value - min_anomaly_value) / 1000.0 + Double.MIN_VALUE;
            int[] detectedTime = new int[24];
            int[] falseAlarmTime = new int[24];
            
            int detected_num, false_alarm_num;
            while (threshold <= max_anomaly_value){
                detected_num = 0; false_alarm_num = 0;
                detectedTime = new int[24];
                falseAlarmTime = new int[24];
                for (EvaluateResult evalRes : routeResult){
                	int evalResTime = Integer.parseInt(evalRes.seg_id.split("_")[1]);
                    if (evalRes.acc_label >= 1 && evalRes.anomaly_val >= threshold){
                    	//System.out.println("anomaly_val: "+evalRes.anomaly_val+" Threshold: "+threshold);
                        detected_num++;
                        detectedTime[evalResTime]++;
                    }
                    if (evalRes.acc_label == 0 && evalRes.anomaly_val >= threshold){
                        false_alarm_num++;
                        falseAlarmTime[evalResTime]++;
                    }
                }
                //System.out.println("Detected: "+detected_num+" Total: "+acc_num);
                //System.out.println("False Alarmed: "+false_alarm_num+" Total: "+total_num);
                double DR = (double) detected_num / (double) acc_num;
                double FAR = (double) false_alarm_num / (double) (total_num);
                
                fw.write(FAR+"\t"+DR+"\t");
                for(int i=0;i<24;i++){
                	double DRt = (double)detectedTime[i]/(double)accTime[i];
                	double FARt = (double) falseAlarmTime[i]/(double)totalTime[i];
                	fw.write(FARt+"\t"+DRt+"\t");
                	System.out.println(i);
                }
                fw.write("\n");
                threshold += step;
                //step *= 2;
            }
            
            System.out.println("4 loop finished");*/
            fw = new FileWriter("data/incident_detection/TW_"+TIME_WINDOW+"/detect_result_Weighted_k="+k+".txt",true);
            //fw.write("Weighted KL Result\n");
            int[] detectedTime = new int[24];
            int[] falseAlarmTime = new int[24];
            //fw.write("FAR\tDR\n");
            double threshold = min_anomaly_value_weighted;
            double step = (max_anomaly_value_weighted - min_anomaly_value_weighted) / 1000.0;
            int detected_num,false_alarm_num, ttd;
            int stepFlag = 0;
            while (threshold <= max_anomaly_value_weighted){
                detected_num = 0; false_alarm_num = 0;
                detectedTime = new int[24];
                falseAlarmTime = new int[24];
                ttd = 0;
                int i = 0;
                for (EvaluateResult evalRes : routeResult){
                	int evalResTime = Integer.parseInt(evalRes.seg_id.split("_")[1]);
                	if (evalRes.isStart == 1){
                		int j = i;
                		while(routeResult.get(j).acc_label == 1 && routeResult.get(j).seg_id.split("_")[0].equals(evalRes.seg_id.split("_")[0])){
                			if (routeResult.get(j).anomaly_val_weighted >= threshold){
                				ttd += (j - i) * 5;
                            	detected_num++;
                                detectedTime[evalResTime]++;
                                break;
                            }
                			j++;
                		}	
                	}
        
                    if (evalRes.acc_label == 0 && evalRes.anomaly_val_weighted >= threshold){
                        false_alarm_num++;
                        falseAlarmTime[evalResTime]++;
                    }
                    i++;
                }
                double DR = (double) detected_num / (double) acc_num;
                double FAR = (double) false_alarm_num / (double) (normal_num);
                double MTTD = (double) ttd/ (double) detected_num;
                fw.write(FAR+"\t"+DR+"\t" + MTTD + "\t");
                System.out.println(k+"\t" +FAR+"\t"+DR+"\t" + MTTD);
                for(int j=0;j<24;j++){
                	double DRt = (double)detectedTime[j]/(double)accTime[j];
                	double FARt = (double) falseAlarmTime[j]/(double)totalTime[j];
                	fw.write(FARt+"\t"+DRt+"\t");
                }
                //System.out.println("Threshold: "+threshold);
                fw.write("\n");
                threshold += step;
                if(FAR < 0.3 && stepFlag == 0) {
                	step += step;
                	stepFlag++;
                }
                if(FAR < 0.2 && stepFlag == 1) {
                	step += 5 * step;
                	stepFlag++;
                }
                if(FAR < 0.15 && stepFlag == 2) {
                	step += step;
                	stepFlag++;
                }
                if(FAR < 0.1 && stepFlag == 3) {
                	step += step;
                	stepFlag++;
                }
                if(FAR < 0.05 && stepFlag == 4) {
                	step += step;
                	stepFlag++;
                }
            }
            System.out.println("all loop finished");
            
            fw.close();
            br.close();
            is.close();
        //}
    }
    public static void main(String[] args) throws FileNotFoundException, IOException{
    	//sun.arch.data.model=64;
    	System.out.println(System.getProperty("sun.arch.data.model")) ;
    	System.out.println("Execution Start");
        /*for (int i=1; i<=4; i++){
        	System.out.println("Fold "+ Integer.toString(i));
            IncidentDetection id = new IncidentDetection();
            id.parseDistribution("data/estimate_result/fold_"+i+"/estimate.out");
            id.parseLabeledRecords("data/label_data/fold_"+i+"/label_result", i);
        }*/
    	int[] klist = {2, 3, 4, 5, 8};
    	
    	for (int t=1; t<9 ; t++){
    		TIME_WINDOW = t;
    		for (int k: klist){
    			System.out.println("Time Window: "+t+" K= "+k);
	    		IncidentDetection id = new IncidentDetection();
	    		id.parseDistribution("data/estimate_result/km/estimate_k"+k+".out");
	    		id.parseLabeledRecords("data/label_data/Labeled_12m_km.txt", k);
    		}
    	}
    	System.out.println("Execution Finished");
    }
}
