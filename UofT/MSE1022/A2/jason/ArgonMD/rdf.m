function RDF=rdf(N,t,x,L,RDF)

G=gr(N,x,L);

for i=1:500
    RDF(i,2) = RDF(i,2) + G(i,2);
end



