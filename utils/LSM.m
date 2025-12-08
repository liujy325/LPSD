function [parameters,rmse]=LSM(match1,match2,change_form)
% This function gets the transformation parameters based on two sets of
% matches using the least square method
match1_xy=match1(:,[1 2]);
match2_xy=match2(:,[1 2]); %generate the initial point set randomly

A=[];
for i=1:size(match1_xy,1)
    A=[A;repmat(match1_xy(i,:),2,2)];
end
B=[1 1 0 0;0 0 1 1];B=repmat(B,size(match1_xy,1),1);
A=A.*B;
B=[1 0;0 1];B=repmat(B,size(match1_xy,1),1);
A=[A B];
t_match2_xy=match2_xy';b=t_match2_xy(:);

if(strcmp(change_form,'affine'))
    [Q,R]=qr(A);
    parameters=R\(Q'*b);
    parameters(7:8)=0;
    N=size(match1,1);
    match1_test=match1(:,[1 2]);match2_test=match2(:,[1 2]);
    M=[parameters(1) parameters(2);parameters(3) parameters(4)];
    match1_test_trans=M*match1_test'+repmat([parameters(5);parameters(6)],1,N);
    match1_test_trans=match1_test_trans';
    test=match1_test_trans-match2_test;
    rmse=sqrt(sum(sum(test.^2))/N);
elseif(strcmp(change_form,'perspective'))
    temp_1=[];
    for i=1:size(match1_xy,1)
    temp_1=[temp_1;repmat(match1_xy(i,:),2,1)];
    end
    temp_1=-1.*temp_1;
    temp_2=repmat(b,1,2);
    temp=temp_1.*temp_2;
    A=[A temp];
    [Q,R]=qr(A);
    parameters=R\(Q'*b);
    N=size(match1,1);
    match1_test=match1_xy';
    match1_test=[match1_test;ones(1,N)];
    M=[parameters(1),parameters(2),parameters(5);...
        parameters(3),parameters(4),parameters(6);...
        parameters(7),parameters(8),1];
    match1_test_trans=M*match1_test;
    match1_test_trans_12=match1_test_trans(1:2,:);
    match1_test_trans_3=match1_test_trans(3,:);
    match1_test_trans_3=repmat(match1_test_trans_3,2,1);
    match1_test_trans=match1_test_trans_12./match1_test_trans_3;
    match1_test_trans=match1_test_trans';
    match2_test=match2_xy;
    test=match1_test_trans-match2_test;
    rmse=sqrt(sum(sum(test.^2))/N);
elseif(strcmp(change_form,'similarity'))
    A=[];
    for i=1:1:size(match1_xy,1)
        A=[A;match1_xy(i,:),1,0;match1_xy(i,2),-match1_xy(i,1),0,1];
    end
    [Q,R]=qr(A);
    parameters=R\(Q'*b);
    parameters(7:8)=0;
    parameters(5:6)=parameters(3:4);
    parameters(3)=-parameters(2);
    parameters(4)=parameters(1);
    N=size(match1,1);
    match1_test=match1(:,[1 2]);match2_test=match2(:,[1 2]);
    M=[parameters(1) parameters(2);parameters(3) parameters(4)];
    match1_test_trans=M*match1_test'+repmat([parameters(5);parameters(6)],1,N);
    match1_test_trans=match1_test_trans';
    test=match1_test_trans-match2_test;
    rmse=sqrt(sum(sum(test.^2))/N);
end                                      
end
