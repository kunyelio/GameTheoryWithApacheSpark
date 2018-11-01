package allocation;

import org.apache.spark.ml.linalg.DenseMatrix;
import org.apache.spark.ml.linalg.DenseVector;
import org.apache.spark.ml.linalg.Vector;
import org.apache.spark.ml.linalg.Vectors;
import java.util.ArrayList;

public class Algorithm {
    public static void main(String[] args){
        Example2_Step2();
    }

    private static void Example1(){
        double cons_set1_vals[] = {0,	8,	4,	2,	8,	4,	2,	8,	4};
        DenseVector consSet1 = new DenseVector(cons_set1_vals);

        double u1_vals[] = {0,	2,	3,	5,	2,	3,	6,	3,	6};
        DenseVector u1 = new DenseVector(u1_vals);


        double cons_set2_vals[] = {4,	4,	0,	0,	8,	6,	0,	4,	2};
        DenseVector consSet2 = new DenseVector(cons_set2_vals);

        double u2_vals[] = {2,	3,	0,	0,	2,	4,	0,	3,	5};
        DenseVector u2 = new DenseVector(u2_vals);


        double cons_set3_vals[] = {1,	2,	4,	2,	4,	2,	1,	8,	4};
        DenseVector consSet3 = new DenseVector(cons_set3_vals);

        double u3_vals[] = {1,	2,	3,	3,	1,	2,	3,	1,	1};
        DenseVector u3 = new DenseVector(u3_vals);


        double cons_set4_vals[] = {2,	1,	0,	0,	2,	1,	0,	1,	0};
        DenseVector consSet4 = new DenseVector(cons_set4_vals);

        double u4_vals[] = {5,	2,	0,	0,	3,	2,	0,	4,	0};
        DenseVector u4 = new DenseVector(u4_vals);

        ArrayList<Vector> Ul = new ArrayList<Vector>();
        Ul.add(u1);
        Ul.add(u2);
        Ul.add(u3);
        Ul.add(u4);

        ArrayList<Vector> consumptionSets = new ArrayList<Vector>();
        consumptionSets.add(consSet1);
        consumptionSets.add(consSet2);
        consumptionSets.add(consSet3);
        consumptionSets.add(consSet4);

        double S_vals[] = {4,	2,	1,	1,	4,	4,	1,	4,	2};
        DenseVector S = new DenseVector(S_vals);

        AssignmentsAndPrice Optimal = getOptimalAssignmentsAndPrice(Ul, S, consumptionSets);
        getAllocations(Optimal.price, Ul, S, consumptionSets, Optimal.assignmentsL);
    }

    private static void Example2_Step2(){

        double cons_set1_vals[] = {0,	0,	100,	50,	0,	0,	0,	0,	10,	10,	30};
        DenseVector consSet1 = new DenseVector(cons_set1_vals);

        double u1_vals[] = {0,	0,	4,	4,	0,	0,	0,	0,	4,	4,	4};
        DenseVector u1 = new DenseVector(u1_vals);


        double cons_set2_vals[] = {0,	0,	0,	0,	0,	0,	0,	0,	40,	40,	20};
        DenseVector consSet2 = new DenseVector(cons_set2_vals);

        double u2_vals[] = {0,	0,	0,	0,	0,	0,	0,	0,	5,	5,	5};
        DenseVector u2 = new DenseVector(u2_vals);


        double cons_set3_vals[] = {0,	0,	0,	0,	0,	0,	0,	0,	20,	20,	10};
        DenseVector consSet3 = new DenseVector(cons_set3_vals);

        double u3_vals[] = {0,	0,	0,	0,	0,	0,	0,	0,	4,	4,	4};
        DenseVector u3 = new DenseVector(u3_vals);


        double cons_set4_vals[] = {50,	50,	0,	0,	100,	50,	50,	100,	30,	30,	0};
        DenseVector consSet4 = new DenseVector(cons_set4_vals);

        DenseVector u4 = (Vectors.zeros(11)).toDense();

        ArrayList<Vector> Ul = new ArrayList<Vector>();
        Ul.add(u1);
        Ul.add(u2);
        Ul.add(u3);
        Ul.add(u4);

        ArrayList<Vector> consumptionSets = new ArrayList<Vector>();
        consumptionSets.add(consSet1);
        consumptionSets.add(consSet2);
        consumptionSets.add(consSet3);
        consumptionSets.add(consSet4);

        double S_vals[] = {50,50,100,50,100,50,50,100,100,100,50};
        DenseVector S = new DenseVector(S_vals);

        AssignmentsAndPrice Optimal = getOptimalAssignmentsAndPrice(Ul, S, consumptionSets);
        getAllocations(Optimal.price, Ul, S, consumptionSets, Optimal.assignmentsL);

    }

    private static void Example2_Step1(){

        double cons_set1_vals[] = {50,	20,	30,	100,	0,	0,	0,	0,	0,	0,	0};
        DenseVector consSet1 = new DenseVector(cons_set1_vals);

        double u1_vals[] = {3,	3,	0,	6,	0,	0,	0,	0,	0,	0,	0};
        DenseVector u1 = new DenseVector(u1_vals);


        double cons_set2_vals[] = {40,	40,	20,	50,	0,	0,	0,	0,	0,	0,	0};
        DenseVector consSet2 = new DenseVector(cons_set2_vals);

        double u2_vals[] = {5,	5,	0,	8,	0,	0,	0,	0,	0,	0,	0};
        DenseVector u2 = new DenseVector(u2_vals);


        double cons_set3_vals[] = {0,	0,	0,	0,	20,	20,	10,	50,	0,	0,	0};
        DenseVector consSet3 = new DenseVector(cons_set3_vals);

        double u3_vals[] = {0,	0,	0,	0,	5,	5,	5,	5,	0,	0,	0};
        DenseVector u3 = new DenseVector(u3_vals);


        double cons_set4_vals[] = {10,10,0,0,100,50,50,100,	30,30,0};
        DenseVector consSet4 = new DenseVector(cons_set4_vals);

        DenseVector u4 = (Vectors.zeros(11)).toDense();

        ArrayList<Vector> Ul = new ArrayList<Vector>();
        Ul.add(u1);
        Ul.add(u2);
        Ul.add(u3);
        Ul.add(u4);

        ArrayList<Vector> consumptionSets = new ArrayList<Vector>();
        consumptionSets.add(consSet1);
        consumptionSets.add(consSet2);
        consumptionSets.add(consSet3);
        consumptionSets.add(consSet4);

        double S_vals[] = {50,	50,	100, 50, 100, 50, 50, 100, 100, 100, 50};
        DenseVector S = new DenseVector(S_vals);

        AssignmentsAndPrice Optimal = getOptimalAssignmentsAndPrice(Ul, S, consumptionSets);
        getAllocations(Optimal.price, Ul, S, consumptionSets, Optimal.assignmentsL);

    }

    private static AssignmentsAndPrice getOptimalAssignmentsAndPrice(ArrayList<Vector> Ul, Vector S, ArrayList<Vector> consumptionSets){
        Vector P = Vectors.zeros(S.size());

        Vector P_new = null;
        AssignmentsAndPrice Optimal = null;
        String[] Omega = base2(P.size());
        while(true){
            AssignmentsAndPrice result = getCMin(Ul,P,consumptionSets,S,Omega);
            P_new = result.price;

            if(P.equals(P_new)){
                Optimal = new AssignmentsAndPrice();
                Optimal.price = P_new;
                Optimal.assignmentsL = result.assignmentsL;
                break;
            }else{
                P = P_new;
            }
        }
        return Optimal;
    }

    private static AssignmentsAndPrice getCMin(ArrayList<Vector> Ul, Vector P,
                                                        ArrayList<Vector> consumptionSets, Vector S,
                                                        String[] Omega){
        ArrayList<double[]> optimal = null;
        double minimum = -1d;
        Vector P_tmp = null;
        String optimalRepresentation = null;

        for(String s:Omega){
            AssignmentsAndUtility tmpObj = null;
            double[] tmp = new double[s.length()];
            for(int j = 0; j < tmp.length; j++){
                tmp[j] = Character.getNumericValue(s.charAt(j));
            }
            DenseVector omega_v = new DenseVector(tmp);
            Vector P_omega = addOrSubtract(P,omega_v,true);
            double maximum = ((new DenseMatrix(1, P_omega.size(), P_omega.toArray())).multiply(S).toArray())[0];
            ArrayList<double[]> tmpArr = new ArrayList<double[]>();
            for(int i = 0; i < Ul.size(); i++){
                tmpObj = getMaxUtility(Ul.get(i), P_omega, consumptionSets.get(i));
                tmpArr.add(tmpObj.assignments);
                maximum +=  tmpObj.utility;
            }

            if(minimum < 0 || maximum < minimum){
                minimum = maximum;
                optimal = tmpArr;
                optimalRepresentation = s;
                P_tmp = P_omega;
            }

        }

        AssignmentsAndPrice result = new AssignmentsAndPrice();
        result.price = P_tmp;
        result.assignmentsL = optimal;
        return result;
    }

    private static ArrayList<DenseVector> getAllocations(Vector P, ArrayList<Vector> Ul, Vector S, ArrayList<Vector>
            consumptionSets, ArrayList<double[]> assignmentsL){
        int K = S.size(); // size of supply vector
        int N = assignmentsL.size(); // size of bidders

        double P_vals[] = P.toArray();

        ArrayList<double[]> Ul_vals = new ArrayList<double[]>();
        for(Vector v:Ul){
            Ul_vals.add(v.toArray().clone());
        }

        ArrayList<double[]> consumptionSet_vals = new ArrayList<double[]>();
        for(Vector v:consumptionSets){
            consumptionSet_vals.add(v.toArray().clone());
        }

        double S_vals[] = S.toArray();

        ArrayList<double[]> Alloc = new ArrayList<double[]>();


        // Initialize allocations
        for(int i = 0; i < N; i++){
            double[] vals = new double[K];
            Alloc.add(vals);
        }

        // For each supply item, perform distribution according to allocations and consumption sets
        for(int j = 0; j < K; j++){
            double Z_j = S_vals[j];
            while(Z_j > 0) {
                int bidderIndex = 0;

                // if for that i there is any bidder with nonzero allocation and ui > price allocate to it
                for (double[] bidder : Alloc) {
                    if (assignmentsL.get(bidderIndex)[j] > 0 && consumptionSet_vals.get(bidderIndex)[j] > 0) {
                        if(Z_j <= 0){
                            break;
                        }
                        if(Ul_vals.get(bidderIndex)[j] - P_vals[j] > 0){
                            double Delta = Math.min(Z_j,consumptionSet_vals.get(bidderIndex)[j]);
                            consumptionSet_vals.get(bidderIndex)[j] = consumptionSet_vals.get(bidderIndex)[j] - Delta;
                            Z_j = Z_j - Delta;
                            bidder[j] = bidder[j] + Delta;
                        }
                    }
                    bidderIndex++;
                }

                //

                bidderIndex = 0;

                for (double[] bidder : Alloc) {
                    if (assignmentsL.get(bidderIndex)[j] > 0 && consumptionSet_vals.get(bidderIndex)[j] > 0) {
                        if(Z_j <= 0){
                            break;
                        }
                        double Delta = Math.min(Z_j,consumptionSet_vals.get(bidderIndex)[j]);
                        consumptionSet_vals.get(bidderIndex)[j] = consumptionSet_vals.get(bidderIndex)[j] - Delta;
                        Z_j = Z_j - Delta;
                        bidder[j] = bidder[j] + Delta;
                    }
                    bidderIndex++;
                }
            }
        }

        ArrayList<DenseVector> allocations = new ArrayList<DenseVector>();
        for(double[] allocation:Alloc){
            allocations.add(new DenseVector(allocation));
        }
        return allocations;
    }

    private static AssignmentsAndUtility getMaxUtility(Vector U, Vector P, Vector consumptionSet){
        double[] tmp1 = addOrSubtract(U,P,false).toArray();
        double[] tmp2 = new double[tmp1.length];
        double[] normalized = new double[tmp1.length];

        int j = 0; // j is an index to resource types
        for(double dbl:tmp1){
            if(dbl < 0d){
                tmp2[j] = 0d;
                normalized[j] = 0d;
                j++;
            }else{
                tmp2[j] = dbl;
                normalized[j] = 1d;
                j++;
            }
        }
        DenseMatrix uMinusp = new DenseMatrix(1, tmp2.length, tmp2);
        double result = (uMinusp.multiply(consumptionSet).toArray())[0];

        AssignmentsAndUtility returnObject = new AssignmentsAndUtility();
        returnObject.utility = result;
        returnObject.assignments = normalized;
        return returnObject;
    }

    private static DenseVector addOrSubtract(Vector a, Vector b, boolean isAdd){
        double[] tmp1 = a.toArray();
        double[] tmp2 = b.toArray();
        double[] tmp3 = new double[tmp1.length];
        int factor = 1;
        if(!isAdd){
            factor = -1;
        }
        for(int i = 0; i < tmp1.length; i++){
            tmp3[i] = tmp1[i] + factor*tmp2[i];
        }
        return new DenseVector(tmp3);
    }

    private static String[] base2(int size){
        int dimension = (int)Math.pow(2,size);
        String[] retArr = new String[dimension];
        for(int i = 0; i < dimension; i++){
            String tmp = Integer.toString(i,2);
            retArr[i] = pad(size,tmp);
        }
        return retArr;
    }


    private static String pad(int size, String suffix){
        String zero = "0";
        StringBuilder bldr = new StringBuilder();
        for(int i = 0; i < size - suffix.length(); i++){
            bldr.append(zero);
        }
        return bldr.toString() + suffix;
    }
}

class AssignmentsAndUtility {
    double[] assignments;
    double utility;
}

class AssignmentsAndPrice{
    ArrayList<double[]> assignmentsL;
    Vector price;
}
