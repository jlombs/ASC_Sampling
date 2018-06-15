function ai = myNChooseKm(n,k)
%%Function which returns the binomial coefficient, assuming double inputs
bi = n - k;
ai = bi + 1;
for i = 2:k
    ai = ai + (ai * bi) / i;
end
end