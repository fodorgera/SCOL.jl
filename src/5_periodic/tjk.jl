function tjk(j, k, τ, Coll)
    tjk_eta = η -> τ / (Coll.div * (Coll.xs[end] - Coll.xs[1])) * η + (j - 1) * τ + k * τ / Coll.div - τ * Coll.xs[end] / (Coll.div * (Coll.xs[end] - Coll.xs[1]))
    # this is valid for divs
    return tjk_eta
end