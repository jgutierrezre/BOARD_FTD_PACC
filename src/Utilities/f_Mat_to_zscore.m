function [ tf_zscore ] = f_Mat_to_zscore( tf_Data )
%F_MAT_TO_ZSCORE normalisa la matriz de tiempo frecuencia con zscore
%   tf_zscore: matriz de tiempo frecuencia normalizada zscore
%   tf_Data: matriz de tiempo frecuencia datos en filas

mean_tf = mean(tf_Data')';
std_tf = std(tf_Data,[],2);
tf_zs = tf_Data-repmat(mean_tf,1,size(tf_Data,2));
tf_zscore = tf_zs./repmat(std_tf,1,size(tf_Data,2));
end

