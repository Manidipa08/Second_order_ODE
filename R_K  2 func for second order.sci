function y=RKsecond(x1,y1,x,f,h)
    y(:,1)=y1
    n=length(x)
    disp("length : ",n)
    for i=1:n-1
        k1=h*f(x(i),y(:,i))
        k2=h*f((x(i)+h),(y(:,i)+k1))
        y(:,i+1)=y(:,i)+(1/2)*(k1+k2)
    end
endfunction

