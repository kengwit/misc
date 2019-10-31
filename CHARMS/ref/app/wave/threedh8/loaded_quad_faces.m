function [fconn]=loaded_quad_faces(fen, hexconn, box, tol)
fen=[fen zeros(size(fen,1), 1)];
for i=1:size(fen, 1)
    fen(i,5)=fen_in_box(box, tol, fen(i,2:4));
end

fconn=[];
for i=1:size(hexconn, 1)
    fconn=[fconn; quad_face(fen(hexconn(i,:),1),fen(hexconn(i,:),5))];
end


function result = fen_in_box(box, tol, xyz)
result=true;
for j=1:3
    result = result & (xyz(j) >= box(j)-tol & xyz(j) <= box(j+3)+tol);
end
        
function f=quad_face(fenid, in_box)
fs=[ 1     4     3     2;...
     1     2     6     5;...
     2     3     7     6;...
     3     4     8     7;...
     4     1     5     8;...
     5     6     7     8];
f=[];
for i=1:6
    if (sum(in_box(fs(i,:))) == 4)
        f=[f; fenid(fs(i,:))'];
    end
end

