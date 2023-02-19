function y=RKfourth(x1,y1,x,f,h)
    y(:,1)=y1
    n= length(x)
    disp("length : ",n)
    for i=1:n-1
        k1=h*f(x(i),y(:,i))
        k2=h*f((x(i)+(h/2)),(y(:,i)+(k1/2)))
        k3=h*f((x(i)+(h/2)),(y(:,i)+(k2/2)))
        k4=h*f((x(i)+h),(y(:,i)+k3))
        y(:,i+1)=y(:,i)+((1/6)*(k1+(2*k2)+(2*k3)+k4))
    end
endfunction

