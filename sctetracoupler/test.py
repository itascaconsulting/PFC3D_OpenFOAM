import sctetracoupler
reload(sctetracoupler)


model = sctetracoupler.ScTetraCoupler("ini-wall")

model.execute()