package bgt.IncidentDetection.Models;

/**
 * Created by admin on 2016/12/2.
 */
public class LabeledRecord {
    public double time;
    public int speed;
    public int acc_flag; // 1 for accident, 0 for no accident
    public String seg_id;
    public int day;
    public int isStart;
    public LabeledRecord (String seg_id, int day, double time, int speed, int acc_flag){
        this.seg_id = seg_id;
        this.day = day;
    	this.time = time;
        this.speed = speed;
        this.acc_flag = acc_flag;
    }
}
