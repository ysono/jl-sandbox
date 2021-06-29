### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 15febb00-0dd8-4272-ae40-16ed6e9847d4
begin
	import Pkg

#     Pkg.activate(mktempdir())
#     Pkg.add([
# 			Pkg.PackageSpec(name="Plots"),
#           Pkg.PackageSpec(name="PlutoUI"),
# 			Pkg.PackageSpec(name="ViscousFlow"),
# 			])

	Pkg.activate(pwd())
	Pkg.instantiate()
end

# ╔═╡ 8fa7dc2f-40f5-4ff8-be41-ce9569f66e95
using Plots, PlutoUI, ViscousFlow

# ╔═╡ f107fe7e-43a7-4a88-860f-6c9d78532108
md"### Problem specification"

# ╔═╡ 3630b8fd-2739-4ebc-b234-75bec2f4d50b
begin
	Re = 200
	U = 1.0
	U∞ = (U, 0.0)
end

# ╔═╡ c746eb6a-248a-48c9-8e4c-2cd382c2fbe5
begin
	xlim = (-4.5, 1.5)
	ylim = (-1.0, 0.5)
	Δx, Δt = setstepsizes(Re, gridRe=4.0)
	Δx, Δt
end

# ╔═╡ 63aa3251-dc40-4690-a62d-bb501579c4a7
md"### Set up bodies"

# ╔═╡ b4126d9a-e5d9-49c4-b3a7-0982c0342e81
function make_teardrop(diameter::Real, l1::Real, θ_slope::Real, gap::Real)
	radius = diameter / 2
	θθ_arc = -π/2 : gap/radius : π/2
	xx_arc = cos.(θθ_arc) * radius
	yy_arc = sin.(θθ_arc) * radius
	
	xx_btm = xx_arc[1]-l1 : gap : xx_arc[1]
	yy_btm = [yy_arc[1] for _ in xx_btm]
	
	x_topleft = xx_btm[1] + diameter / tan(θ_slope)

	xx_slope = x_topleft : -gap * cos(θ_slope) : xx_btm[1]
	yy_slope = LinRange(yy_arc[end], yy_btm[1], size(xx_slope, 1))
	
	xx_top = xx_arc[end] : -gap : xx_slope[1]
	yy_top = [radius for _ in xx_top]
	
	xx = vcat(xx_top, xx_slope, xx_btm, xx_arc)
	yy = vcat(yy_top, yy_slope, yy_btm, yy_arc)
	
	BasicBody(xx, yy)
end

# ╔═╡ 29c6f875-851d-4b2b-bc6b-0df7d5336b04
plot(make_teardrop(3, 9, π/6, 0.3), size=(300, 200))

# ╔═╡ 6dab7a8d-d8f7-4759-8523-ef08579ab01b
function make_edge(diameter::Real, θ_slope::Real, repeat::Int, gap::Real)
	function pattern()
		radius = diameter / 2
		θθ_arc = -π/2 : gap/radius : π/2+θ_slope
		xx_arc = cos.(θθ_arc) * radius
		yy_arc = sin.(θθ_arc) * radius
		
		x_left = xx_arc[end] - (yy_arc[end] - yy_arc[1]) / tan(θ_slope)
		
		xx_itv = xx_arc[end] : -gap * cos(θ_slope) : x_left
		yy_itv = yy_arc[end] : -gap * sin(θ_slope) : yy_arc[1]
		
		vcat(xx_arc, xx_itv), vcat(yy_arc, yy_itv)
	end
	
	pad = diameter / 2 * 1.2
	
	xx, yy = pattern()
	for _ in 2 : max(repeat, 1)
		xx_, yy_ = pattern()
		append!(xx, xx_ .+ xx[end] .- pad)
		append!(yy, yy_)
	end
	
	append!(xx, [xx[end], pad, pad])
	append!(yy, [pad, pad, yy[1]])
	
	xx_last = xx[end] : -gap : xx[1]
	yy_last = [yy[1] for _ in xx_last]
	append!(xx, xx_last)
	append!(yy, yy_last)
	
	BasicBody(xx, yy)
end

# ╔═╡ ebecf77c-034a-4783-a2c0-820f9397a643
plot(make_edge(5, π/6, 2, 0.3))

# ╔═╡ 9d9b9556-7717-46cd-a94a-2aa178e81079
function make_scene()
	bl = BodyList()
	tl = RigidTransformList()
	
	rpt = 2
	
	function edge(flip::Bool=false)
		e = make_edge(0.5, π/12, rpt, 1.5Δt)
		if flip; e.y .*= -1; end
		e
	end
	function tear(flip::Bool=false)
		t = make_teardrop(0.3, .8, π/6, 1.5Δt)
		if flip; t.y .*= -1; end
		t
	end
	
	push!(bl, edge())
	push!(tl, RigidTransform((0., 0.), 0.))
	
	push!(bl, edge(true))
	push!(tl, RigidTransform((1.1, -0.45), 0.))
	
	for offset in 0 : rpt-1
	# for offset in 0 : 1
	# for offset in 0 : 0
		push!(bl, tear(true))
		push!(tl, RigidTransform((float(offset) * -2.2, 0.), π/12))
		
		push!(bl, tear())
		push!(tl, RigidTransform((float(offset) * -2.2 + 1.1, -0.45), -π/12))
	end
	
	tl(bl)
	
	bl
	
	# tear(true)
end

# ╔═╡ 1cce7d87-b239-4f86-9936-fd38ceb7ba06
begin
	body = make_scene()
	plot(body, xlim=xlim, ylim=ylim)
	# plot(body)
end

# ╔═╡ b48b1f59-8514-46a6-8fe2-f610ac58806c
md"### Construct the system structure"

# ╔═╡ 6619de6b-7156-449d-bf27-ef9f8476d507
sys=NavierStokes(Re,Δx,xlim,ylim,Δt,body,freestream=U∞)

# ╔═╡ 19d0125f-e3ae-4baf-bbc0-767b2ec8518e
md"### Initialize"

# ╔═╡ ec86b599-02f4-45a8-ae99-2538fde7de27
u0 = newstate(sys)

# ╔═╡ a47a0991-1b31-44ea-8c04-9f6e4bb7370a
begin
	tspan = (0.0, 3.0)
	integrator = init(u0,tspan,sys)
end

# ╔═╡ 95135b6a-876f-4265-b742-db3bcb6ed9d3
md"### Solve"

# ╔═╡ 92faff30-5b7e-403b-88c6-f5abd45bdb48
step!(integrator,3.0)

# ╔═╡ d7678688-7b46-449a-8785-d6f2d5f0abe9
plot(
	plot(vorticity(integrator), sys, levels=100, title="Vorticity"),
	plot(streamfunction(integrator), sys, title="Streamlines"),
	plot(pressure(integrator), sys, title="Pressure"),
	plot(velocity(integrator), sys, levels=100, title="Velocity"),
	layout=(4,1),
	size=(800, 1200)
)

# ╔═╡ Cell order:
# ╠═15febb00-0dd8-4272-ae40-16ed6e9847d4
# ╠═8fa7dc2f-40f5-4ff8-be41-ce9569f66e95
# ╟─f107fe7e-43a7-4a88-860f-6c9d78532108
# ╠═3630b8fd-2739-4ebc-b234-75bec2f4d50b
# ╠═c746eb6a-248a-48c9-8e4c-2cd382c2fbe5
# ╟─63aa3251-dc40-4690-a62d-bb501579c4a7
# ╠═b4126d9a-e5d9-49c4-b3a7-0982c0342e81
# ╠═29c6f875-851d-4b2b-bc6b-0df7d5336b04
# ╠═6dab7a8d-d8f7-4759-8523-ef08579ab01b
# ╠═ebecf77c-034a-4783-a2c0-820f9397a643
# ╠═9d9b9556-7717-46cd-a94a-2aa178e81079
# ╠═1cce7d87-b239-4f86-9936-fd38ceb7ba06
# ╟─b48b1f59-8514-46a6-8fe2-f610ac58806c
# ╠═6619de6b-7156-449d-bf27-ef9f8476d507
# ╟─19d0125f-e3ae-4baf-bbc0-767b2ec8518e
# ╠═ec86b599-02f4-45a8-ae99-2538fde7de27
# ╠═a47a0991-1b31-44ea-8c04-9f6e4bb7370a
# ╟─95135b6a-876f-4265-b742-db3bcb6ed9d3
# ╠═92faff30-5b7e-403b-88c6-f5abd45bdb48
# ╠═d7678688-7b46-449a-8785-d6f2d5f0abe9
