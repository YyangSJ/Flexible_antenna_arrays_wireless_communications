function AE = CosineAnt(theta,phi,ka)

if cos(phi)<0
    AE=0;
else
    AE=sqrt(2*(ka+1)*(sin(theta)*cos(phi))^ka);
end

end

