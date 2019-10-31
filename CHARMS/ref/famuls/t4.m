child_map=[ 0 4 6 7;
            4 1 5 8;
            6 5 2 9;
            7 8 9 3;
            8 7 9 4;
            8 5 4 9;
            6 4 5 9;
            9 6 4 7 ];
child_map=child_map+1; % shift by 1
child_map_reflected=[  0 4 6 7;
  4 1 5 8;
  6 5 2 9;
  7 8 9 3;
  4 7 8 6;
  4 5 6 8;
  9 6 5 8;
  9 8 7 6 ];
child_map_reflected= child_map_reflected+1;

rst=[   0    0    0  ;
        1    0    0  ;
        0    1    0  ;
        0    0    1  ;
       1/2   0    0  ;
       1/2  1/2   0  ;
        0   1/2   0  ;
        0    0   1/2 ;
       1/2   0   1/2 ;
        0   1/2  1/2 ];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the child to parent map
%
fprintf (1, '// Map child to parent\n');
fprintf (1, '// ##########################################\n');
cm = child_map;
for child=1:8

   A=zeros(12,12);
   rhs=zeros(12);

   for i=1:4
      A((i-1)*3+1,1:3)=rst(i,:);
      A((i-1)*3+2,4:6)=rst(i,:);
      A((i-1)*3+3,7:9)=rst(i,:);
      A((i-1)*3+1:(i-1)*3+3,10:12)=eye(3,3);
      rhs((i-1)*3+1:(i-1)*3+3)=rst(cm(child,i),:);
   end
   x=A\rhs;
   B=zeros(3,3);
   S=zeros(3);
   B(1,:)=x(1:3);
   B(2,:)=x(4:6);
   B(3,:)=x(7:9);
   S=x(10:12);
%   disp('B matrix'); B;
 %  disp('Shift '); S;
   child_rst=[rst(1,:);
              rst(2,:);
              rst(3,:);
              rst(4,:)];
   parent_rst=B*child_rst'+[S;S;S;S]';
   Binv=inv(B);
   child_rst=Binv*(parent_rst-[S;S;S;S]');
   fprintf (1, '{{{%g,%g,%g},{%g,%g,%g},{%g,%g,%g}},{%g,%g,%g}}, // child %d\n',B(1,1),B(1,2),B(1,3),B(2,1),B(2,2),B(2,3),B(3,1),B(3,2),B(3,3),S(1),S(2),S(3),child-1);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the child to parent map (reflected)
%
fprintf (1, '// Map child to parent (reflected)\n');
fprintf (1, '// ##########################################\n');
cm = child_map_reflected;
for child=1:8

   A=zeros(12,12);
   rhs=zeros(12);

   for i=1:4
      A((i-1)*3+1,1:3)=rst(i,:);
      A((i-1)*3+2,4:6)=rst(i,:);
      A((i-1)*3+3,7:9)=rst(i,:);
      A((i-1)*3+1:(i-1)*3+3,10:12)=eye(3,3);
      rhs((i-1)*3+1:(i-1)*3+3)=rst(cm(child,i),:);
   end
   x=A\rhs;
   B=zeros(3,3);
   S=zeros(3);
   B(1,:)=x(1:3);
   B(2,:)=x(4:6);
   B(3,:)=x(7:9);
   S=x(10:12);
%   disp('B matrix'); B;
 %  disp('Shift '); S;
   child_rst=[rst(1,:);
              rst(2,:);
              rst(3,:);
              rst(4,:)];
   parent_rst=B*child_rst'+[S;S;S;S]';
   Binv=inv(B);
   child_rst=Binv*(parent_rst-[S;S;S;S]');
   fprintf (1, '{{{%g,%g,%g},{%g,%g,%g},{%g,%g,%g}},{%g,%g,%g}}, // child %d\n',B(1,1),B(1,2),B(1,3),B(2,1),B(2,2),B(2,3),B(3,1),B(3,2),B(3,3),S(1),S(2),S(3),child-1);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the  parent to child map
%
fprintf (1, '// Map parent to child\n');
fprintf (1, '// ##########################################\n');
cm = child_map;
for child=1:8

   A=zeros(12,12);
   rhs=zeros(12);

   for i=1:4
      A((i-1)*3+1,1:3)=rst(i,:);
      A((i-1)*3+2,4:6)=rst(i,:);
      A((i-1)*3+3,7:9)=rst(i,:);
      A((i-1)*3+1:(i-1)*3+3,10:12)=eye(3,3);
      rhs((i-1)*3+1:(i-1)*3+3)=rst(cm(child,i),:);
   end
   x=A\rhs;
   B=zeros(3,3);
   S=zeros(3);
   B(1,:)=x(1:3);
   B(2,:)=x(4:6);
   B(3,:)=x(7:9);
   S=x(10:12);
%   disp('B matrix'); B;
 %  disp('Shift '); S;
   parent_rst=[rst(cm(child,1),:);
              rst(cm(child,2),:);
              rst(cm(child,3),:);
              rst(cm(child,4),:)];
   child_rst=B*parent_rst'+[S;S;S;S]';
   Binv=inv(B);
   parent_rst=Binv*child_rst-[S;S;S;S]';
   Sinv=-(Binv*S');
   fprintf (1, '{{{%g,%g,%g},{%g,%g,%g},{%g,%g,%g}},{%g,%g,%g}}, // child %d\n',Binv(1,1),Binv(1,2),Binv(1,3),Binv(2,1),Binv(2,2),Binv(2,3),Binv(3,1),Binv(3,2),Binv(3,3),Sinv(1),Sinv(2),Sinv(3),child-1);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the  parent to child map
%
fprintf (1, '// Map parent to child (reflected)\n');
fprintf (1, '// ##########################################\n');
cm = child_map_reflected;
for child=1:8

   A=zeros(12,12);
   rhs=zeros(12);

   for i=1:4
      A((i-1)*3+1,1:3)=rst(i,:);
      A((i-1)*3+2,4:6)=rst(i,:);
      A((i-1)*3+3,7:9)=rst(i,:);
      A((i-1)*3+1:(i-1)*3+3,10:12)=eye(3,3);
      rhs((i-1)*3+1:(i-1)*3+3)=rst(cm(child,i),:);
   end
   x=A\rhs;
   B=zeros(3,3);
   S=zeros(3);
   B(1,:)=x(1:3);
   B(2,:)=x(4:6);
   B(3,:)=x(7:9);
   S=x(10:12);
%   disp('B matrix'); B;
 %  disp('Shift '); S;
   parent_rst=[rst(cm(child,1),:);
              rst(cm(child,2),:);
              rst(cm(child,3),:);
              rst(cm(child,4),:)];
   child_rst=B*parent_rst'+[S;S;S;S]';
   Binv=inv(B);
   parent_rst=Binv*child_rst-[S;S;S;S]';
   Sinv=-(Binv*S');
   fprintf (1, '{{{%g,%g,%g},{%g,%g,%g},{%g,%g,%g}},{%g,%g,%g}}, // child %d\n',Binv(1,1),Binv(1,2),Binv(1,3),Binv(2,1),Binv(2,2),Binv(2,3),Binv(3,1),Binv(3,2),Binv(3,3),Sinv(1),Sinv(2),Sinv(3),child-1);

end
